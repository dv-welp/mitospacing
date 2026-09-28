import os
import re
import sys
from datetime import datetime

import numpy as np
import tifffile
import roifile
from scipy.ndimage import map_coordinates
from scipy.signal import find_peaks

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox, CheckButtons

from openpyxl import Workbook, load_workbook

import tkinter as tk
from tkinter import messagebox

# ---------------- CONFIG ----------------
ROOT = r"path/to/your/date_folder""
PX_PER_UM = 6.59
PEAK_DISTANCE = 5
PEAK_PROMINENCE = 10
CLICK_TOLERANCE_UM = 5
PEAK_LABEL_PREFIX = "mit#"               
IMAGE_EXTENSIONS = ('.tif', '.tiff')
ROI_SUFFIX_HINTS = ('roi',)              # words to strip when matching roi filenames
IMAGE_SUFFIX_HINTS = ('mosaic',)         # words to strip when matching image filenames (also catches "MosaicJ", "mosaic2", etc.)
LOG_SHEET_NAME = "corrections_log"
UNMATCHED_SHEET_NAME = "unmatched_files"
UTILITY_SHEETS = {LOG_SHEET_NAME, UNMATCHED_SHEET_NAME}   # excluded when analysis looks for genotype sheets
# -------------------------------------------------------

DATE_NAME = os.path.basename(ROOT.rstrip("\\/")) or "data"
WORKBOOK_PATH = os.path.join(ROOT, f"mitochondrial_distribution_{DATE_NAME}.xlsx")
CLICK_TOLERANCE_PX = CLICK_TOLERANCE_UM * PX_PER_UM

# Guarantee backspace works as a "back to previous view" shortcut (same action as the
# toolbar's back arrow), regardless of what this matplotlib version defaults to.
# Only triggers when focus isn't in a text field (e.g. the Comment box), so backspace
# still deletes text normally while you're typing a comment.
if 'backspace' not in plt.rcParams['keymap.back']:
    plt.rcParams['keymap.back'].append('backspace')


# ==================== profile extraction ====================

def get_profile_and_length(tif_path, roi_path):
    img = tifffile.imread(tif_path)
    roi = roifile.ImagejRoi.fromfile(roi_path)
    coords = roi.coordinates()
    xs, ys = coords[:, 0], coords[:, 1]

    profile = []
    coords_xy = []  # (x, y) image coordinate for every profile sample, same order/index as profile
    for i in range(len(xs) - 1):
        x1, y1, x2, y2 = xs[i], ys[i], xs[i + 1], ys[i + 1]
        seg_len = int(np.hypot(x2 - x1, y2 - y1))
        if seg_len == 0:
            continue
        seg_x = np.linspace(x1, x2, seg_len)
        seg_y = np.linspace(y1, y2, seg_len)
        profile.extend(map_coordinates(img, [seg_y, seg_x], order=1))
        coords_xy.extend(zip(seg_x, seg_y))

    profile = np.array(profile)
    length_um = len(profile) / PX_PER_UM
    return profile, length_um, coords_xy, img


def _nearest_coord(coords_xy, index_float):
    """Rounds a (possibly fractional) profile index to the nearest sample and
    returns its stored (x, y) image coordinate."""
    idx = int(round(index_float))
    idx = max(0, min(idx, len(coords_xy) - 1))
    return coords_xy[idx]


def _local_perpendicular(coords_xy, index_float, delta=5):
    """Returns (point_xy, unit_perpendicular) at a profile index -- the
    direction 90 degrees off the ROI's local tangent there, used to offset
    peak-indicator arrows so they never sit on top of the peak itself."""
    idx = int(round(index_float))
    idx = max(0, min(idx, len(coords_xy) - 1))
    i0 = max(0, idx - delta)
    i1 = min(len(coords_xy) - 1, idx + delta)
    p0 = np.array(coords_xy[i0])
    p1 = np.array(coords_xy[i1])
    point = np.array(coords_xy[idx])

    tangent = p1 - p0
    norm = np.linalg.norm(tangent)
    if norm == 0:
        tangent_unit = np.array([1.0, 0.0])
    else:
        tangent_unit = tangent / norm
    perp = np.array([-tangent_unit[1], tangent_unit[0]])
    return point, perp


# ==================== folder discovery ====================

def normalize_name(name):
    """Lowercase, strip known image/roi hint words (plus anything glued directly
    onto them, e.g. 'mosaic' matches 'MosaicJ', 'mosaic2', etc.), collapse whitespace.
    Matching only -- never affects what's written to the workbook."""
    base, _ext = os.path.splitext(name)
    base = base.lower()
    for hint in ROI_SUFFIX_HINTS + IMAGE_SUFFIX_HINTS:
        base = re.sub(re.escape(hint) + r"\w*", "", base)
    base = re.sub(r"\s+", " ", base).strip()
    return base


def discover_animals(genotype_path):
    """
    Returns (animals, unmatched).
      animals: list of {"image_id", "tif_path", "roi_path"}
      unmatched: list of {"filename", "reason"} -- flagged for manual review

    Handles both layouts automatically:
      - nested: one subfolder per animal, each with exactly one image + one .roi
      - flat: images and .roi files sitting directly in genotype_path,
              matched by normalized filename
    """
    entries = os.listdir(genotype_path)
    subdirs = [e for e in entries if os.path.isdir(os.path.join(genotype_path, e))]
    files = [e for e in entries if os.path.isfile(os.path.join(genotype_path, e))]

    animals = []
    unmatched = []

    # --- nested layout ---
    for sub in sorted(subdirs):
        sub_path = os.path.join(genotype_path, sub)
        try:
            tif_name = next(f for f in os.listdir(sub_path) if f.lower().endswith(IMAGE_EXTENSIONS))
            roi_name = next(f for f in os.listdir(sub_path) if f.lower().endswith('.roi'))
            animals.append({
                "image_id": sub,
                "tif_path": os.path.join(sub_path, tif_name),
                "roi_path": os.path.join(sub_path, roi_name),
            })
        except StopIteration:
            unmatched.append({"filename": sub, "reason": "subfolder missing an image or a .roi file"})

    # --- flat layout ---
    tif_files = [f for f in files if f.lower().endswith(IMAGE_EXTENSIONS)]
    roi_files = [f for f in files if f.lower().endswith('.roi')]

    roi_by_key = {}
    for r in roi_files:
        roi_by_key.setdefault(normalize_name(r), []).append(r)

    used_rois = set()
    for t in sorted(tif_files):
        key = normalize_name(t)
        candidates = roi_by_key.get(key, [])
        if len(candidates) == 1:
            roi_name = candidates[0]
            used_rois.add(roi_name)
            animals.append({
                "image_id": key,
                "tif_path": os.path.join(genotype_path, t),
                "roi_path": os.path.join(genotype_path, roi_name),
            })
        elif len(candidates) == 0:
            unmatched.append({"filename": t, "reason": "no matching .roi found (exact normalized match)"})
        else:
            unmatched.append({"filename": t, "reason": f"multiple possible .roi matches: {candidates}"})

    for r in roi_files:
        if r not in used_rois:
            unmatched.append({"filename": r, "reason": "no matching image found (exact normalized match)"})

    return animals, unmatched


# ==================== workbook helpers ====================

def init_workbook():
    if os.path.exists(WORKBOOK_PATH):
        return load_workbook(WORKBOOK_PATH)
    wb = Workbook()
    wb.remove(wb.active)  # sheets are created per genotype as needed
    log_ws = wb.create_sheet(LOG_SHEET_NAME)
    log_ws.append(["timestamp", "genotype", "image", "action", "position_px", "position_um"])
    unmatched_ws = wb.create_sheet(UNMATCHED_SHEET_NAME)
    unmatched_ws.append(["genotype", "filename", "reason"])
    return wb


def get_genotype_sheet(wb, genotype_name):
    if genotype_name in wb.sheetnames:
        return wb[genotype_name]
    ws = wb.create_sheet(genotype_name)
    ws.append(["Id"])
    return ws


def find_image_column(ws, image_id):
    for col in range(2, ws.max_column + 1):
        if ws.cell(row=1, column=col).value == image_id:
            return col
    return None


def read_existing_state(ws, col):
    if col is None:
        return [], None, False
    comment = None
    reviewed = False
    peaks_um = []
    for row in range(2, ws.max_row + 1):
        label = ws.cell(row=row, column=1).value
        val = ws.cell(row=row, column=col).value
        if label == "comment":
            comment = val
        elif label == "reviewed":
            reviewed = bool(val)
        elif isinstance(label, str) and label.startswith(PEAK_LABEL_PREFIX) and isinstance(val, (int, float)):
            peaks_um.append(val)
    return sorted(peaks_um), comment, reviewed


def write_image_column(ws, image_id, peaks_um, length_um, comment, reviewed):
    col = find_image_column(ws, image_id)
    if col is None:
        col = ws.max_column + 1
        ws.cell(row=1, column=col, value=image_id)

    labels_needed = ["length_um", "density", "number", "comment", "reviewed"]
    label_to_row = {}
    for row in range(2, ws.max_row + 1):
        val = ws.cell(row=row, column=1).value
        if val:
            label_to_row[val] = row

    next_row = ws.max_row + 1
    for label in labels_needed:
        if label not in label_to_row:
            ws.cell(row=next_row, column=1, value=label)
            label_to_row[label] = next_row
            next_row += 1

    n_peaks = len(peaks_um)
    density = (n_peaks / length_um * 100) if length_um else None

    ws.cell(row=label_to_row["length_um"], column=col, value=round(length_um, 2))
    ws.cell(row=label_to_row["density"], column=col, value=round(density, 3) if density else None)
    ws.cell(row=label_to_row["number"], column=col, value=n_peaks)
    ws.cell(row=label_to_row["comment"], column=col, value=comment or "")
    ws.cell(row=label_to_row["reviewed"], column=col, value=bool(reviewed))

    for row in range(2, ws.max_row + 1):
        val = ws.cell(row=row, column=1).value
        if isinstance(val, str) and val.startswith(PEAK_LABEL_PREFIX):
            ws.cell(row=row, column=col, value=None)

    for i, p in enumerate(sorted(peaks_um), start=1):
        peak_label = f"{PEAK_LABEL_PREFIX}{i}"
        if peak_label not in label_to_row:
            row = ws.max_row + 1
            ws.cell(row=row, column=1, value=peak_label)
            label_to_row[peak_label] = row
        ws.cell(row=label_to_row[peak_label], column=col, value=round(p, 2))


def log_correction(wb, genotype_name, image_id, action, position_px):
    log_ws = wb[LOG_SHEET_NAME]
    log_ws.append([
        datetime.now().isoformat(timespec='seconds'),
        genotype_name,
        image_id,
        action,
        round(position_px, 1),
        round(position_px / PX_PER_UM, 2),
    ])


def log_unmatched(wb, genotype_name, unmatched_list):
    if not unmatched_list:
        return
    ws = wb[UNMATCHED_SHEET_NAME]
    existing = {(row[0].value, row[1].value) for row in ws.iter_rows(min_row=2) if row[0].value}
    for item in unmatched_list:
        key = (genotype_name, item["filename"])
        if key not in existing:
            ws.append([genotype_name, item["filename"], item["reason"]])


# ==================== annotation panel ====================

def annotate_animal(wb, genotype_name, animal):
    image_id = animal["image_id"]
    profile, length_um, coords_xy, img = get_profile_and_length(animal["tif_path"], animal["roi_path"])

    ws = get_genotype_sheet(wb, genotype_name)
    col = find_image_column(ws, image_id)
    existing_peaks_um, existing_comment, _ = read_existing_state(ws, col)

    if existing_peaks_um:
        peaks_px = [p * PX_PER_UM for p in existing_peaks_um]
        print(f"[{genotype_name}] {image_id}: resuming from {len(peaks_px)} previously saved peaks")
    else:
        raw_peaks, _ = find_peaks(profile, distance=PEAK_DISTANCE, prominence=PEAK_PROMINENCE)
        peaks_px = raw_peaks.tolist()
        print(f"[{genotype_name}] {image_id}: starting from {len(peaks_px)} automatically detected peaks")

    session_corrections = []

    fig, (ax_img, ax_profile) = plt.subplots(1, 2, figsize=(15, 6))
    plt.subplots_adjust(left=0.03, right=0.98, top=0.86, bottom=0.22, wspace=0.1)

    # ---- image panel ----
    ax_img.imshow(img, cmap='gray')
    ax_img.set_xticks([])
    ax_img.set_yticks([])

    roi_xs = [c[0] for c in coords_xy]
    roi_ys = [c[1] for c in coords_xy]
    roi_line, = ax_img.plot(roi_xs, roi_ys, color='cyan', linewidth=1, alpha=0.8)

    hover_marker, = ax_img.plot([], [], marker='o', markerfacecolor='none', markeredgecolor='yellow',
                                markeredgewidth=2, markersize=9, zorder=10)
    hover_marker.set_visible(False)

    ARROW_OFFSET = 10   # px away from the peak the arrow starts
    ARROW_GAP = 3        # px short of the peak the arrowhead stops
    peak_arrows = []      # current arrow annotation artists, redrawn on every change

    def redraw_peak_arrows():
        for artist in peak_arrows:
            artist.remove()
        peak_arrows.clear()
        for p in peaks_px:
            point, perp = _local_perpendicular(coords_xy, p)
            start = point + perp * ARROW_OFFSET
            end = point + perp * ARROW_GAP
            arrow = ax_img.annotate(
                "", xy=tuple(end), xytext=tuple(start),
                arrowprops=dict(arrowstyle="->", color="lime", lw=1.5),
            )
            peak_arrows.append(arrow)

    redraw_peak_arrows()

    # ---- ROI show/hide toggle ----
    roi_check_ax = plt.axes([0.02, 0.89, 0.13, 0.05], frameon=False)
    roi_checkbox = CheckButtons(roi_check_ax, ["Show ROI"], [True])

    def on_roi_toggle(label):
        roi_line.set_visible(not roi_line.get_visible())
        fig.canvas.draw_idle()

    roi_checkbox.on_clicked(on_roi_toggle)

    # ---- profile panel ----
    ax_profile.plot(profile)
    scatter = ax_profile.scatter(
        peaks_px,
        [profile[int(round(p))] for p in peaks_px] if peaks_px else [],
        facecolors='none', edgecolors='green', linewidths=1.5, zorder=5,
    )
    ax_profile.set_xlabel("position along ROI (px)")

    # ---- header: title and length/count info kept on separate rows, clear gap ----
    fig.suptitle(f"{genotype_name} / {image_id}", fontsize=13, fontweight='bold', y=0.97)
    info_text = fig.text(0.5, 0.90, "", ha='center', fontsize=10)

    def update_info_text():
        info_text.set_text(f"length: {length_um:.1f} um    |    mitochondria: {len(peaks_px)}")

    update_info_text()

    textbox_ax = plt.axes([0.15, 0.04, 0.4, 0.06])
    textbox = TextBox(textbox_ax, "Comment: ", initial=existing_comment or "")

    def redraw():
        if peaks_px:
            ys = [profile[int(round(p))] for p in peaks_px]
            scatter.set_offsets(np.c_[peaks_px, ys])
        else:
            scatter.set_offsets(np.empty((0, 2)))
        redraw_peak_arrows()
        update_info_text()
        fig.canvas.draw_idle()

    def on_click(event):
        if event.inaxes != ax_profile or event.xdata is None:
            return
        x_click = event.xdata
        if peaks_px:
            diffs = [abs(x_click - p) for p in peaks_px]
            min_diff = min(diffs)
            min_idx = diffs.index(min_diff)
        else:
            min_diff = None

        if min_diff is not None and min_diff <= CLICK_TOLERANCE_PX:
            removed = peaks_px.pop(min_idx)
            session_corrections.append(("remove", removed))
        else:
            peaks_px.append(x_click)
            peaks_px.sort()
            session_corrections.append(("add", x_click))
        redraw()

    def on_motion(event):
        if event.inaxes == ax_profile and event.xdata is not None:
            x_img, y_img = _nearest_coord(coords_xy, event.xdata)
            hover_marker.set_data([x_img], [y_img])
            hover_marker.set_visible(True)
        else:
            hover_marker.set_visible(False)
        fig.canvas.draw_idle()

    def save_current_state(mark_reviewed):
        peaks_um_final = [p / PX_PER_UM for p in peaks_px]
        comment_text = textbox.text
        write_image_column(ws, image_id, peaks_um_final, length_um, comment_text, mark_reviewed)
        for action, pos_px in session_corrections:
            log_correction(wb, genotype_name, image_id, action, pos_px)
        session_corrections.clear()

        try:
            wb.save(WORKBOOK_PATH)
            status = "DONE" if mark_reviewed else "in progress"
            print(f"Saved [{genotype_name}] {image_id} ({len(peaks_um_final)} peaks, {status}).")
        except PermissionError:
            fallback_path = WORKBOOK_PATH.replace(
                ".xlsx", f"_UNSAVED_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"
            )
            print("\n" + "!" * 70)
            print(f"COULD NOT SAVE -- {WORKBOOK_PATH} is open elsewhere (likely in Excel).")
            print(f"Your edits for {image_id} were NOT lost -- saved instead to:")
            print(f"  {fallback_path}")
            print("Close the workbook everywhere else, then manually merge this backup")
            print("in, or just rerun the script once the file is free.")
            print("!" * 70 + "\n")
            try:
                wb.save(fallback_path)
            except Exception as e:
                print(f"Backup save ALSO failed ({e}) -- do not close this window until you've "
                      f"copied down your peak positions manually.")

    state = {"mark_reviewed": False}

    def on_save_button(event):
        state["mark_reviewed"] = False
        plt.close(fig)

    def on_done_button(event):
        state["mark_reviewed"] = True
        plt.close(fig)

    def on_close(event):
        save_current_state(state["mark_reviewed"])

    save_ax = plt.axes([0.6, 0.04, 0.15, 0.06])
    save_button = Button(save_ax, "Save & Close")
    save_button.on_clicked(on_save_button)

    done_ax = plt.axes([0.78, 0.04, 0.17, 0.06])
    done_button = Button(done_ax, "Mark Done & Close")
    done_button.on_clicked(on_done_button)

    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('motion_notify_event', on_motion)
    fig.canvas.mpl_connect('close_event', on_close)

    plt.show()


# ==================== analysis ====================

def parse_and_verify_sheet(filepath, sheet_name):
    """
    Reads one genotype sheet from the workbook produced by this pipeline.
    Finds 'density', 'length_um', 'number' by label (order-independent), and
    treats every row whose label starts with PEAK_LABEL_PREFIX ("mit#") as a
    coordinate row -- so it doesn't matter how many peaks a given image has.
    """
    import pandas as pd

    df = pd.read_excel(filepath, sheet_name=sheet_name, header=None)
    labels = df[0].astype(str).str.strip().str.lower().tolist()

    try:
        density_row_idx = labels.index('density')
        length_row_idx = labels.index('length_um')
        number_row_idx = labels.index('number')
    except ValueError:
        missing = [label for label in ['density', 'length_um', 'number'] if label not in labels]
        raise ValueError(f"Sheet '{sheet_name}' is missing these expected row labels in column A: {missing}")

    id_row_idx = 0  # "Id" is always the header row in this workbook's format

    sample_data = df.iloc[density_row_idx, 1:]
    n_samples = len(sample_data.dropna())

    sample_ids = df.iloc[id_row_idx, 1:n_samples + 1].tolist()

    excel_densities = pd.to_numeric(df.iloc[density_row_idx, 1:n_samples + 1], errors='coerce').values
    process_lengths = pd.to_numeric(df.iloc[length_row_idx, 1:n_samples + 1], errors='coerce').values
    excel_counts = pd.to_numeric(df.iloc[number_row_idx, 1:n_samples + 1], errors='coerce').values

    for idx, sample_id in enumerate(sample_ids):
        if pd.isna(process_lengths[idx]) or process_lengths[idx] == 0:
            continue
        calc_density = (excel_counts[idx] / process_lengths[idx]) * 100
        if not np.isclose(calc_density, excel_densities[idx], atol=1e-2, equal_nan=True):
            print(f"[MISMATCH WARNING] Sheet: {sheet_name} | Sample: {sample_id} | "
                  f"Excel Density: {excel_densities[idx]:.3f}, Calculated: {calc_density:.3f}")

    peak_row_indices = [i for i, label in enumerate(labels)
                        if isinstance(label, str) and label.startswith(PEAK_LABEL_PREFIX)]
    coord_matrix = df.iloc[peak_row_indices, 1:n_samples + 1].apply(pd.to_numeric, errors='coerce').values

    raw_absolute_distances = []
    for col_idx in range(n_samples):
        coords = coord_matrix[:, col_idx]
        coords = coords[~np.isnan(coords) & (coords > 0)]
        raw_absolute_distances.append(list(sorted(coords)))

    return sample_ids, excel_densities, process_lengths, raw_absolute_distances


def compute_spatial_metrics(raw_absolute_distances, comp_sizes_10, comp_sizes_20, no_of_parcel=5):
    """
    Computes IMD tables, binned distributions, Fano factors,
    and Morisita's Index for both 10um and 20um steps.
    """
    intermitochondrial_distances = []
    first_mit_distances = []

    for coords in raw_absolute_distances:
        if len(coords) == 0:
            intermitochondrial_distances.append([])
            first_mit_distances.append(np.nan)
            continue
        first_mit_distances.append(coords[0])
        imds = [coords[0]]
        if len(coords) > 1:
            imds.extend(np.diff(coords))
        intermitochondrial_distances.append(imds)

    bins = np.arange(0, 66 + 3, 3)
    abs_freqs, norm_freqs, fano_factors = [], [], []

    for imds in intermitochondrial_distances:
        counts, _ = np.histogram(imds, bins=bins)
        abs_freqs.append(counts)
        norm_freqs.append(counts / len(imds) if len(imds) > 0 else counts * 0)

        mean_c = np.mean(counts)
        var_c = np.var(counts, ddof=1)
        fano_factors.append(var_c / mean_c if mean_c > 0 else np.nan)

    def run_morisita_loop(sizes):
        results = np.zeros((len(sizes), len(raw_absolute_distances)))
        for idx, comp_size in enumerate(sizes):
            size_of_parcel = comp_size / no_of_parcel
            for j, coords in enumerate(raw_absolute_distances):
                parcel_counts = np.zeros(no_of_parcel)
                for n in range(no_of_parcel):
                    parcel_counts[n] = np.sum((np.array(coords) > size_of_parcel * n) &
                                               (np.array(coords) <= size_of_parcel * (n + 1)))
                N = np.sum(parcel_counts)
                if N > 1:
                    results[idx, j] = no_of_parcel * (np.sum(parcel_counts * (parcel_counts - 1)) / (N * (N - 1)))
                else:
                    results[idx, j] = np.nan
        return results

    morisita_10 = run_morisita_loop(comp_sizes_10)
    morisita_20 = run_morisita_loop(comp_sizes_20)

    return (intermitochondrial_distances, first_mit_distances,
            np.array(abs_freqs), np.array(norm_freqs), fano_factors, morisita_10, morisita_20)


def save_table_figure(df, title, out_path):
    """Renders a DataFrame as a simple table image and saves it."""
    n_rows, n_cols = len(df), len(df.columns)
    fig_width = max(6, 1.6 * n_cols)
    fig_height = 0.5 * n_rows + 1.2
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis('off')
    table = ax.table(
        cellText=df.round(4).astype(str).values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center',
        bbox=[0, 0, 1, 1],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.auto_set_column_width(col=list(range(n_cols)))
    ax.set_title(title, fontweight='bold', pad=18, fontsize=14)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()


def generate_publication_plots(results_store, genotypes, comp_range_10, comp_range_20, output_dir):
    import seaborn as sns

    def outpath(name):
        return os.path.join(output_dir, name)

    palette = sns.color_palette("Set2", len(genotypes))

    all_densities = [results_store[g]['densities'] for g in genotypes]
    all_lengths = [results_store[g]['lengths'] for g in genotypes]
    all_first_mit = [results_store[g]['first_mit'] for g in genotypes]

    morphology_traits = [
        (all_densities, "Mitochondrial Density", r"Density (per 100 $\mu$m)", "plot_density.png"),
        (all_lengths, "Process Length", r"Length ($\mu$m)", "plot_length.png"),
        (all_first_mit, "Distance of 1st Mit from CB", r"Distance ($\mu$m)", "plot_first_mit.png")
    ]

    for data_list, title, ylabel, filename in morphology_traits:
        plt.figure(figsize=(7, 6))

        sns.violinplot(data=data_list, palette=palette, inner=None, alpha=0.3, linewidth=1.5)
        sns.boxplot(data=data_list, palette=palette, width=0.12, boxprops=dict(alpha=0.7), showfliers=False, color="black")
        sns.stripplot(data=data_list, color='black', alpha=0.2, jitter=0.12, size=4)

        for i, d in enumerate(data_list):
            clean_d = np.array(d)[~np.isnan(d)]
            mean_val = np.mean(clean_d)
            sem_val = np.std(clean_d, ddof=1) / np.sqrt(len(clean_d)) if len(clean_d) > 0 else 0

            plt.errorbar(i, mean_val, yerr=sem_val, fmt='o', color='darkred', elinewidth=2.5, markersize=8, zorder=5)

            y_axis_range = plt.gca().get_ylim()[1] - plt.gca().get_ylim()[0]
            y_text_pos = plt.gca().get_ylim()[1] - (y_axis_range * 0.08)

            plt.text(i, y_text_pos, f"{mean_val:.2f} \u00b1 {sem_val:.2f}",
                     ha='center', va='center', color='darkred', weight='bold', fontsize=11,
                     bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', boxstyle='round,pad=0.2'))

        plt.xticks(range(len(genotypes)), genotypes)
        plt.ylabel(ylabel, labelpad=10)
        plt.title(title, fontweight='bold', pad=20, fontsize=14)
        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        plt.savefig(outpath(filename), dpi=300)
        plt.close()

    bins_3um = np.arange(0, 66 + 3, 3)
    bin_labels_upper = [str(bins_3um[i + 1]) for i in range(len(bins_3um) - 1)]
    bin_centers = (bins_3um[:-1] + bins_3um[1:]) / 2

    for i, g in enumerate(genotypes):
        plt.figure(figsize=(8, 5))
        nf = results_store[g]['norm_freq']

        mean_profile = np.mean(nf, axis=0)
        sem_profile = np.std(nf, axis=0, ddof=1) / np.sqrt(len(nf))

        plt.bar(bin_centers, mean_profile, width=2.4, yerr=sem_profile,
                color=palette[i], alpha=0.7, edgecolor=palette[i], linewidth=1.5,
                error_kw=dict(ecolor='black', elinewidth=1.5, capsize=2))

        plt.title(f"Normalised Frequency Histogram - {g}", fontweight='bold', pad=15, fontsize=14)
        plt.xlabel(r"Intermitochondrial Distance Bins ($\mu$m)", labelpad=10)
        plt.ylabel("Probability Density", labelpad=10)

        display_labels = [label if idx % 2 == 1 else "" for idx, label in enumerate(bin_labels_upper)]
        plt.xticks(bin_centers, display_labels, rotation=0, fontsize=10)

        sns.despine(offset=10, trim=True)
        plt.tight_layout()
        plt.savefig(outpath(f"plot_histogram_{g}.png"), dpi=300)
        plt.close()

    plt.figure(figsize=(7, 6))
    for i, g in enumerate(genotypes):
        m10 = results_store[g]['morisita_10']
        plt.errorbar(comp_range_10, np.nanmean(m10, axis=1),
                     yerr=np.nanstd(m10, axis=1, ddof=1) / np.sqrt(m10.shape[1]),
                     fmt='-o', color=palette[i], label=g, linewidth=2, capsize=3, markersize=5)

    plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5)
    plt.xlabel(r'Compartment size ($\mu$m)', labelpad=10)
    plt.ylabel(r'I-delta index ($I_\delta$)', labelpad=10)
    plt.title(r"Morisita Index ($10\,\mu$m steps)", fontweight='bold', pad=15, fontsize=14)
    plt.legend(frameon=True, fontsize=10, loc='best')
    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    plt.savefig(outpath('plot_morisita_10um.png'), dpi=300)
    plt.close()

    plt.figure(figsize=(7, 6))
    for i, g in enumerate(genotypes):
        m20 = results_store[g]['morisita_20']
        plt.errorbar(comp_range_20, np.nanmean(m20, axis=1),
                     yerr=np.nanstd(m20, axis=1, ddof=1) / np.sqrt(m20.shape[1]),
                     fmt='-s', color=palette[i], label=g, linewidth=2, capsize=3, markersize=5)

    plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5)
    plt.xlabel(r'Compartment size ($\mu$m)', labelpad=10)
    plt.ylabel(r'I-delta index ($I_\delta$)', labelpad=10)
    plt.title(r"Morisita Index ($20\,\mu$m steps)", fontweight='bold', pad=15, fontsize=14)
    plt.legend(frameon=True, fontsize=10, loc='best')
    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    plt.savefig(outpath('plot_morisita_20um.png'), dpi=300)
    plt.close()

    print("\n[Success] All plots generated as clean standalone figures!")


def run_analysis(workbook_path):
    """
    Runs the full statistics + plotting pipeline on the given workbook.
    Heavy analysis-only libraries are imported here, not at module level,
    so a plain review session never loads them.
    """
    import pandas as pd
    import seaborn as sns
    from scipy.stats import shapiro, kruskal, f_oneway
    import scikit_posthocs as sp

    sns.set_theme(style="ticks", context="talk")
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['font.family'] = 'sans-serif'

    output_dir = os.path.join(os.path.dirname(workbook_path), "analysis_output")
    os.makedirs(output_dir, exist_ok=True)

    def outpath(name):
        return os.path.join(output_dir, name)

    comp_range_10 = np.arange(80, 410, 10)
    comp_range_20 = np.arange(80, 420, 20)

    xls = pd.ExcelFile(workbook_path)
    genotypes = [s for s in xls.sheet_names if s not in UTILITY_SHEETS]
    print(f"[Info] Found {len(genotypes)} genotype sheets to analyze: {genotypes}")

    results_store = {}
    fano_summary_list = []

    for genotype in genotypes:
        ids, dens, lengths, coords = parse_and_verify_sheet(workbook_path, genotype)
        metrics = compute_spatial_metrics(coords, comp_range_10, comp_range_20)

        results_store[genotype] = {
            'ids': ids, 'densities': dens, 'lengths': lengths,
            'imds': metrics[0], 'first_mit': metrics[1],
            'abs_freq': metrics[2], 'norm_freq': metrics[3], 'fano': metrics[4],
            'morisita_10': metrics[5], 'morisita_20': metrics[6]
        }

        fano_summary_list.append({
            'Genotype': genotype,
            'Fano Factor (Mean)': np.nanmean(metrics[4]),
            'Fano Factor (SEM)': np.nanstd(metrics[4], ddof=1) / np.sqrt(np.sum(~np.isnan(metrics[4])))
        })

    print("\n================== FANO FACTOR METRIC REPORT ==================")
    fano_df = pd.DataFrame(fano_summary_list)
    print(fano_df.to_string(index=False))

    output_filename = outpath('Genotypes_Troubleshooting_Summary.xlsx')
    with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
        summary_rows = []
        for g in genotypes:
            res = results_store[g]
            summary_rows.append({
                'Genotype': g,
                'Density Mean': np.mean(res['densities']), 'Density SEM': np.std(res['densities'], ddof=1) / np.sqrt(len(res['densities'])),
                'Length Mean': np.mean(res['lengths']), 'Length SEM': np.std(res['lengths'], ddof=1) / np.sqrt(len(res['lengths'])),
                '1st Mit Dist Mean': np.nanmean(res['first_mit']), '1st Mit Dist SEM': np.nanstd(res['first_mit'], ddof=1) / np.sqrt(np.sum(~np.isnan(res['first_mit'])))
            })
            max_len = max(len(x) for x in res['imds'])
            padded = [x + [np.nan] * (max_len - len(x)) for x in res['imds']]
            pd.DataFrame(np.array(padded).T, columns=res['ids']).to_excel(writer, sheet_name=f'{g}_Step_Distances', index=False)
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Global_Morphology_Summary', index=False)
    print(f"\n[Success] Unified validation sheet saved as '{output_filename}'!")

    # ---- statistics ----
    print("\n================== STATISTICAL HYPOTHESIS REPORT ==================")
    metrics_to_test = {
        'Mitochondrial Density': [results_store[g]['densities'] for g in genotypes],
        'Process Length': [results_store[g]['lengths'] for g in genotypes],
        'Distance to 1st Mit': [[val for val in results_store[g]['first_mit'] if not pd.isna(val)] for g in genotypes]
    }

    stats_report_rows = []
    dunn_matrices = {}

    for metric_name, data_groups in metrics_to_test.items():
        print(f"\nTesting Metric: {metric_name}")

        all_normal = True
        for i, g in enumerate(genotypes):
            clean_group_data = [v for v in data_groups[i] if not pd.isna(v)]

            if len(clean_group_data) < 3:
                print(f"  -> Shapiro-Wilk ({g}): Sample size too small (N={len(clean_group_data)}). Defaulting to Non-Normal.")
                all_normal = False
                continue

            stat_sw, p_val_sw = shapiro(clean_group_data)
            is_group_normal = p_val_sw > 0.05
            print(f"  -> Shapiro-Wilk ({g}): p-value = {p_val_sw:.4f} " + ("(Normal)" if is_group_normal else "(NON-NORMAL)"))

            if not is_group_normal:
                all_normal = False

        cleaned_groups = [[v for v in group if not pd.isna(v)] for group in data_groups]

        if all_normal:
            stat, p_val = f_oneway(*cleaned_groups)
            test_name = "One-Way ANOVA"
            print(f"  Parametric One-Way ANOVA Result: F-stat = {stat:.4f}, p-value = {p_val:.4f}")
            if p_val < 0.05:
                print("  Significant difference identified. Pairwise comparisons should be evaluated via Tukey's HSD.")
        else:
            stat, p_val = kruskal(*cleaned_groups)
            test_name = "Kruskal-Wallis"
            print(f"  Non-Parametric Kruskal-Wallis Result: H-stat = {stat:.4f}, p-value = {p_val:.4f}")

            if p_val < 0.05:
                print("  Significant difference identified. Running Dunn's Post-Hoc Test (Holm-corrected):")
                flat_data = []
                flat_groups = []
                for idx, g in enumerate(genotypes):
                    flat_data.extend(cleaned_groups[idx])
                    flat_groups.extend([g] * len(cleaned_groups[idx]))

                df_stat = pd.DataFrame({'Value': flat_data, 'Group': flat_groups})
                dunn_matrix = sp.posthoc_dunn(df_stat, val_col='Value', group_col='Group', p_adjust='holm')
                print(dunn_matrix.to_string())
                dunn_matrices[metric_name] = dunn_matrix
            else:
                print("  No statistically significant baseline difference between genotypes.")

        stats_report_rows.append({
            'Metric': metric_name,
            'Test': test_name,
            'Statistic': round(stat, 4),
            'p-value': round(p_val, 4),
            'Significant (p<0.05)': 'Yes' if p_val < 0.05 else 'No',
        })

    stats_df = pd.DataFrame(stats_report_rows)

    # ---- stats + fano as table figures ----
    save_table_figure(fano_df, "Fano Factor Summary", outpath("table_fano_factor.png"))
    save_table_figure(stats_df, "Statistical Test Summary", outpath("table_stats_summary.png"))
    for metric_name, dunn_matrix in dunn_matrices.items():
        safe_name = metric_name.replace(" ", "_").lower()
        save_table_figure(
            dunn_matrix.reset_index().rename(columns={'index': 'Group'}),
            f"Dunn's Post-Hoc ({metric_name})",
            outpath(f"table_dunn_{safe_name}.png"),
        )

    # ---- plots ----
    generate_publication_plots(results_store, genotypes, comp_range_10, comp_range_20, output_dir)

    print(f"\n[Success] Analysis complete. Outputs saved in: {output_dir}")


def confirm_and_maybe_run_analysis():
    root = tk.Tk()
    root.withdraw()
    proceed = messagebox.askyesno(
        "Review complete",
        f"All images reviewed.\n\nWorkbook:\n{WORKBOOK_PATH}\n\n"
        "Open it now to check if you like, then click Yes to run analysis, "
        "or No to just stop here.",
    )
    root.destroy()
    if proceed:
        run_analysis(WORKBOOK_PATH)
    else:
        print("Skipped analysis for now -- run with --analyze whenever you're ready.")


# ==================== main ====================

def main():
    if "--analyze" in sys.argv:
        if not os.path.exists(WORKBOOK_PATH):
            print(f"No workbook found at {WORKBOOK_PATH} -- nothing to analyze yet.")
            return
        run_analysis(WORKBOOK_PATH)
        return

    wb = init_workbook()

    genotype_dirs = sorted(
        d for d in os.listdir(ROOT)
        if os.path.isdir(os.path.join(ROOT, d))
    )
    print(f"Found {len(genotype_dirs)} genotype folder(s) under {ROOT}")

    n_processed = 0
    n_skipped = 0
    any_unreviewed_remaining = False

    for genotype in genotype_dirs:
        genotype_path = os.path.join(ROOT, genotype)
        animals, unmatched = discover_animals(genotype_path)
        log_unmatched(wb, genotype, unmatched)
        if unmatched:
            print(f"[{genotype}] {len(unmatched)} file(s) flagged in '{UNMATCHED_SHEET_NAME}' sheet -- please check")

        ws = get_genotype_sheet(wb, genotype)

        for animal in animals:
            image_id = animal["image_id"]
            col = find_image_column(ws, image_id)
            _, _, already_reviewed = read_existing_state(ws, col)

            if already_reviewed:
                n_skipped += 1
                continue

            print(f"\n--- [{genotype}] {image_id} ---")
            annotate_animal(wb, genotype, animal)
            n_processed += 1

            col2 = find_image_column(ws, image_id)
            _, _, now_reviewed = read_existing_state(ws, col2)
            if not now_reviewed:
                any_unreviewed_remaining = True

    print(f"\nRan successfully. Processed {n_processed} image(s), skipped {n_skipped} already-done image(s).")
    print(f"Workbook saved at: {WORKBOOK_PATH}")

    if not any_unreviewed_remaining and (n_processed + n_skipped) > 0:
        confirm_and_maybe_run_analysis()


if __name__ == "__main__":
    main()
