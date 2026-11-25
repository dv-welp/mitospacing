"""
Mitochondrial distribution pipeline -- review + analysis, one file.

Folder layout expected:
    ROOT (a date folder)
      |- genotype_1/
      |     (flat layout: image + .roi files directly here, matched by name)
      |     OR (nested layout: one subfolder per animal, each with one image + one .roi)
      |- genotype_2/
      |- ...

Run modes:
    python mitospacing.py            -> review loop, then offer to run analysis
                                        once every image is marked done
    python mitospacing.py --analyze  -> skip the review loop entirely and just
                                        run analysis on the workbook as it is now

Stopping part-way through reviewing, and resuming later:
    The review loop opens one image at a time and there is no "quit" button:
    "Save & Close", "Mark Done & Close" and the window's X all save the image
    you are on and then the script opens the NEXT one. To stop for the day:
      1. On the image you're on, click "Save & Close" (or "Mark Done & Close"
         if you're finished with it). This writes its peaks + comment to the
         workbook.
      2. When the next image's window appears, leave it alone, click the
         terminal and press Ctrl+C. Nothing is lost: everything you saved is
         already in the workbook, and the window you just opened has no edits.
      (If you press Ctrl+C while an image has edits you have NOT saved with one
       of the buttons, those edits are lost.)
    To resume: run the same command again (python mitospacing.py) with the same
    ROOT set below. Images you marked done are skipped; images you only "Save &
    Close"d re-open with the peaks and comment you left; untouched images start
    from the auto-detected peaks. Make sure the workbook isn't open in Excel
    first, or the save will fall back to an _UNSAVED_ backup file.

Review panel controls:
    - click near a green dot to REMOVE that peak
    - click empty space to ADD a peak there
    - type in the Comment box to leave a note
    - "Save & Close": saves progress, image stays open for later review
    - "Mark Done & Close": saves progress AND marks the image reviewed,
      so future runs skip it
    - closing the window (X) behaves like "Save & Close" -- nothing is
      ever lost by closing, on purpose or by accident
    - if the workbook can't be saved (e.g. it's open in Excel), your edits
      are written to a timestamped backup file instead of being lost

To re-open an image you already marked done: open the workbook, find that
image's column, clear/uncheck its "reviewed" cell, save+close the workbook,
then rerun the script.

Analysis (run_analysis) reads the workbook (one sheet per genotype, one column
per animal) and writes stats tables, a summary workbook and figures to an
"analysis_output" folder next to it. Workbook format it understands
(column A = row labels, column B onward = animals):
    Id            <- header row: one animal ID per column
    length(um) / length_um
    density       <- ignored for stats; recomputed from peak count and length
    number        <- ignored for stats; recomputed from peak count
    comment, reviewed, blank rows, "*" marker rows ... are all ignored
    mit #1, mit #2, ... / mit#1, mit#2, ...   <- peak positions in um from the cell body
Row order, label spelling (spaces, "#", "(um)" vs "_um") and formulas-with-cached-values
are all handled. Sheets named corrections_log / unmatched_files are skipped.
It imports pandas/seaborn/scipy.stats/scikit_posthocs itself, lazily, only when
it actually runs -- a plain review session never pays for loading those libraries.

Needs: numpy, pandas, openpyxl, matplotlib, seaborn, scipy, scikit-posthocs,
       tifffile, roifile
"""

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
ROOT = os.path.dirname(os.path.abspath(__file__))   # the folder this script sits in; to use a different folder, replace with r"C:\path\to\folder"
PX_PER_UM = 6.5923
PEAK_DISTANCE = 5
PEAK_PROMINENCE = 10
CLICK_TOLERANCE_UM = 5
PEAK_LABEL_PREFIX = "mit#"               # rename here if you want a different label later
IMAGE_EXTENSIONS = ('.tif', '.tiff')
ROI_SUFFIX_HINTS = ('roi',)              # words to strip when matching roi filenames
IMAGE_SUFFIX_HINTS = ('mosaic',)         # words to strip when matching image filenames (also catches "MosaicJ", "mosaic2", etc.)
LOG_SHEET_NAME = "corrections_log"
UNMATCHED_SHEET_NAME = "unmatched_files"
UTILITY_SHEETS = {LOG_SHEET_NAME, UNMATCHED_SHEET_NAME}   # excluded when analysis looks for genotype sheets
ANALYSIS_DIR_NAME = "analysis_output"                    # analysis writes here; never treated as a genotype folder
NOT_ITALIC = {"wildtype", "wild type", "wt", "n2"}   # names that are not gene names; remove entries to italicise everything
BIN_WIDTH_UM = 3
BIN_MAX_UM = 66
COMP_RANGE_20 = np.arange(80, 420, 20)               # Morisita compartment sizes (um)
N_PARCELS = 5
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

def _load_analysis_libs():
    """Imports the analysis-only libraries on first use and publishes them as
    module-level names (pd, sns, sp, shapiro, ...) for the functions below, so a
    plain review session never pays for loading them."""
    global pd, sns, sp, shapiro, kruskal, f_oneway, tukey_hsd
    import pandas as pd
    import seaborn as sns
    import scikit_posthocs as sp
    from scipy.stats import shapiro, kruskal, f_oneway, tukey_hsd


def gstyle(name):
    """Font style for a genotype name."""
    return "normal" if str(name).strip().lower() in NOT_ITALIC else "italic"


def style_texts(texts, names):
    for t, n in zip(texts, names):
        t.set_fontstyle(gstyle(n))


# ==================== reading the workbook ====================

def _norm_label(x):
    return re.sub(r"[\s_#()]", "", str(x).strip().lower()) if x is not None else ""


def read_sheet(wb, sheet_name):
    ws = wb[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    labels = [_norm_label(r[0]) if r else "" for r in rows]

    def find(prefix):
        for i, l in enumerate(labels):
            if l.startswith(prefix):
                return i
        return None

    id_i, len_i = find("id"), find("length")
    dens_i, num_i = find("density"), find("number")
    if id_i is None or len_i is None:
        raise ValueError(f"Sheet '{sheet_name}': could not find an 'Id' row and a 'length' row in column A.")
    mit_rows = [i for i, l in enumerate(labels) if re.fullmatch(r"mit\d+", l)]
    if not mit_rows:
        raise ValueError(f"Sheet '{sheet_name}': no 'mit #n' rows found.")
    # Peak block = first 'mit #' row to the bottom of the sheet. Rows in it with a blank
    # label still count as peaks (the number row's COUNT() range includes them).
    peak_rows = [i for i in range(mit_rows[0], len(rows)) if labels[i] == "" or i in mit_rows]
    unlabeled = [i + 1 for i in peak_rows if i not in mit_rows
                 and any(isinstance(v, (int, float)) and not isinstance(v, bool) for v in rows[i][1:])]
    if unlabeled:
        print(f"[Info] {sheet_name}: peak values found in unlabeled rows {unlabeled}; included.")

    def num(v):
        return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else np.nan

    ids, lengths, peaks, stored_dens = [], [], [], []
    n_cols = max(len(r) for r in rows)
    for c in range(1, n_cols):
        sid = rows[id_i][c] if c < len(rows[id_i]) else None
        if sid is None or str(sid).strip() == "":
            continue
        L = num(rows[len_i][c]) if c < len(rows[len_i]) else np.nan
        col_peaks = []
        for r in peak_rows:
            v = rows[r][c] if c < len(rows[r]) else None
            if isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0:
                col_peaks.append(float(v))
        ids.append(str(sid).strip())
        lengths.append(L)
        peaks.append(sorted(col_peaks))
        stored_dens.append(num(rows[dens_i][c]) if dens_i is not None and c < len(rows[dens_i]) else np.nan)

        if num_i is not None:
            stored_n = num(rows[num_i][c]) if c < len(rows[num_i]) else np.nan
            if not np.isnan(stored_n) and int(stored_n) != len(col_peaks):
                print(f"[WARNING] {sheet_name} | {ids[-1]}: 'number' row says {int(stored_n)} "
                      f"but {len(col_peaks)} peak values found.")

    lengths = np.array(lengths)
    counts = np.array([len(p) for p in peaks])
    with np.errstate(divide="ignore", invalid="ignore"):
        density = counts / lengths * 100
    stored_dens = np.array(stored_dens)
    for sid, d, sd in zip(ids, density, stored_dens):
        if not np.isnan(sd) and not np.isclose(d, sd, atol=1e-2):
            print(f"[WARNING] {sheet_name} | {sid}: stored density {sd:.3f} vs recomputed {d:.3f} (using recomputed).")
    return ids, density, lengths, peaks


# ==================== metrics ====================

def compute_spatial_metrics(peaks, comp_sizes, no_of_parcel=N_PARCELS):
    imds_all, first_mit = [], []
    for coords in peaks:
        if len(coords) == 0:
            imds_all.append([])
            first_mit.append(np.nan)
            continue
        first_mit.append(coords[0])
        imds = [coords[0]]
        if len(coords) > 1:
            imds.extend(np.diff(coords))
        imds_all.append(list(imds))

    bins = np.arange(0, BIN_MAX_UM + BIN_WIDTH_UM, BIN_WIDTH_UM)
    abs_freq, norm_freq, fano = [], [], []
    for imds in imds_all:
        counts, _ = np.histogram(imds, bins=bins)
        abs_freq.append(counts)
        norm_freq.append(counts / len(imds) if len(imds) > 0 else counts * 0)
        m = np.mean(counts)
        fano.append(np.var(counts, ddof=1) / m if m > 0 else np.nan)

    morisita = np.full((len(comp_sizes), len(peaks)), np.nan)
    for idx, comp in enumerate(comp_sizes):
        parcel = comp / no_of_parcel
        for j, coords in enumerate(peaks):
            arr = np.array(coords)
            pc = np.array([np.sum((arr > parcel * n) & (arr <= parcel * (n + 1))) for n in range(no_of_parcel)])
            N = pc.sum()
            if N > 1:
                morisita[idx, j] = no_of_parcel * (np.sum(pc * (pc - 1)) / (N * (N - 1)))
    return imds_all, np.array(first_mit), np.array(norm_freq), np.array(fano), morisita


def sem(x):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return np.std(x, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan


# ==================== table figures ====================

def save_table_figure(df, title, out_path, italic_first_col=False, italic_header=False):
    n_rows, n_cols = len(df), len(df.columns)
    fig, ax = plt.subplots(figsize=(max(6, 1.6 * n_cols), 0.5 * n_rows + 1.2))
    ax.axis("off")
    table = ax.table(cellText=df.round(4).astype(str).values, colLabels=df.columns,
                     cellLoc="center", loc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.auto_set_column_width(col=list(range(n_cols)))
    for (r, c), cell in table.get_celld().items():
        if italic_first_col and c == 0 and r > 0:
            cell.get_text().set_fontstyle(gstyle(df.iloc[r - 1, 0]))
        if italic_header and r == 0 and c > 0:
            cell.get_text().set_fontstyle(gstyle(df.columns[c]))
    ax.set_title(title, fontweight="bold", pad=18, fontsize=14)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ==================== plots ====================

def stars(p):
    return "****" if p < 0.0001 else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def sem_note(ax):
    ax.text(1.0, 1.01, "mean \u00b1 SEM", transform=ax.transAxes, ha="right", va="bottom", fontsize=10)


def draw_brackets(ax, genotypes, pmat):
    """Bracket + stars for every significant pair in the post-hoc p-value matrix."""
    if pmat is None:
        return
    pairs = []
    for i in range(len(genotypes)):
        for j in range(i + 1, len(genotypes)):
            p = pmat.loc[genotypes[i], genotypes[j]]
            if stars(p):
                pairs.append((i, j, p))
    if not pairs:
        return
    pairs.sort(key=lambda t: (t[1] - t[0], t[0]))
    lo, hi = ax.get_ylim()
    step = (hi - lo) * 0.10
    tiers = []                                   # list of lists of (i, j) already placed per tier
    for i, j, p in pairs:
        t = 0
        while t < len(tiers) and any(not (j < a or i > b) for a, b in tiers[t]):
            t += 1
        if t == len(tiers):
            tiers.append([])
        tiers[t].append((i, j))
        y = hi + step * (t + 0.5)
        tick = step * 0.18
        ax.plot([i, i, j, j], [y - tick, y, y, y - tick], color="black", linewidth=1.5, clip_on=False)
        ax.text((i + j) / 2, y, stars(p), ha="center", va="bottom", fontsize=13, fontweight="bold")
    ax.set_ylim(lo, hi + step * (len(tiers) + 0.9))


def generate_plots(res, genotypes, output_dir, pvals):
    out = lambda n: os.path.join(output_dir, n)
    pal = dict(zip(genotypes, sns.color_palette("Set2", len(genotypes))))

    # ---- morphology plots: identical figure + axes size ----
    traits = [
        ("densities", "Mitochondrial Density", r"Density (per 100 $\mu$m)", "plot_density.png"),
        ("lengths", "Process Length", r"Length ($\mu$m)", "plot_length.png"),
        ("first_mit", "Distance of First Mit from Cell Body", r"Distance of first mit ($\mu$m)", "plot_first_mit.png"),
    ]
    for key, title, ylabel, fname in traits:
        long = pd.DataFrame([(g, v) for g in genotypes for v in res[g][key] if not np.isnan(v)],
                            columns=["Genotype", "Value"])
        fig = plt.figure(figsize=(7, 6))
        fig.subplots_adjust(left=0.18, right=0.97, bottom=0.14, top=0.88)
        ax = fig.add_subplot(111)
        kw = dict(data=long, x="Genotype", y="Value", order=genotypes, ax=ax)
        sns.violinplot(**kw, hue="Genotype", palette=pal, legend=False, inner=None, linewidth=1.5, alpha=0.3, cut=0)
        sns.boxplot(**kw, hue="Genotype", palette=pal, legend=False, width=0.12, showfliers=False,
                    boxprops=dict(alpha=0.7))
        sns.stripplot(**kw, color="black", alpha=0.2, jitter=0.12, size=4)
        for i, g in enumerate(genotypes):
            vals = np.asarray(res[g][key], dtype=float)
            vals = vals[~np.isnan(vals)]
            ax.errorbar(i, np.mean(vals), yerr=sem(vals), fmt="o", color="darkred", elinewidth=2.5,
                        markersize=8, zorder=5)
        draw_brackets(ax, genotypes, pvals.get(key))
        ax.set_xlabel("")
        ax.set_ylabel(ylabel, labelpad=10)
        ax.set_title(title, fontweight="bold", pad=22, fontsize=14)
        sem_note(ax)
        sns.despine(ax=ax, offset=10, trim=True)
        # tick labels are set LAST: despine(trim=True) rebuilds the ticks and would undo the italics
        ax.set_xticks(range(len(genotypes)))
        ax.set_xticklabels(genotypes)
        for t, g in zip(ax.get_xticklabels(), genotypes):
            t.set_fontstyle(gstyle(g))
        fig.savefig(out(fname), dpi=300)
        plt.close(fig)

    # ---- histograms: one shared y-limit, one fixed axes size ----
    bins = np.arange(0, BIN_MAX_UM + BIN_WIDTH_UM, BIN_WIDTH_UM)
    centers = (bins[:-1] + bins[1:]) / 2
    labels_upper = [str(b) for b in bins[1:]]
    hist = {g: (np.mean(res[g]["norm_freq"], axis=0),
                np.std(res[g]["norm_freq"], axis=0, ddof=1) / np.sqrt(len(res[g]["norm_freq"])))
            for g in genotypes}
    ymax = max(np.nanmax(m + s) for m, s in hist.values()) * 1.08

    for g in genotypes:
        fig = plt.figure(figsize=(8, 5))
        fig.subplots_adjust(left=0.16, right=0.97, bottom=0.22, top=0.82)
        ax = fig.add_subplot(111)
        m, s = hist[g]
        ax.bar(centers, m, width=BIN_WIDTH_UM * 0.8, yerr=s, color=pal[g], alpha=0.7, edgecolor=pal[g],
               linewidth=1.5, error_kw=dict(ecolor="black", elinewidth=1.5, capsize=2))
        fig.suptitle("Normalised Frequency Histogram", fontweight="bold", fontsize=14, y=0.97)
        ax.set_title(g, fontstyle=gstyle(g), fontsize=13, pad=8)
        sem_note(ax)
        ax.set_xlabel(r"Intermitochondrial Distance Bins ($\mu$m)", labelpad=10)
        ax.set_ylabel("Probability Density", labelpad=10)
        ax.set_ylim(0, ymax)
        ax.set_xticks(centers)
        ax.set_xticklabels([l if i % 2 == 1 else "" for i, l in enumerate(labels_upper)], fontsize=10)
        sns.despine(ax=ax, offset=10, trim=True)
        fig.savefig(out(f"plot_histogram_{g}.png"), dpi=300)
        plt.close(fig)

    # ---- Morisita (20 um steps) ----
    fig = plt.figure(figsize=(7, 6))
    fig.subplots_adjust(left=0.18, right=0.96, bottom=0.16, top=0.88)
    ax = fig.add_subplot(111)
    for g in genotypes:
        mor = res[g]["morisita"]
        ax.errorbar(COMP_RANGE_20, np.nanmean(mor, axis=1),
                    yerr=np.nanstd(mor, axis=1, ddof=1) / np.sqrt(np.sum(~np.isnan(mor), axis=1)),
                    fmt="-s", color=pal[g], label=g, linewidth=2, capsize=3, markersize=5)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.5)
    ax.set_xlabel(r"Compartment size ($\mu$m)", labelpad=10)
    ax.set_ylabel(r"I-delta index ($I_\delta$)", labelpad=10)
    ax.set_title(r"Morisita Index ($20\,\mu$m steps)", fontweight="bold", pad=15, fontsize=14)
    sem_note(ax)
    leg = ax.legend(frameon=True, fontsize=10, loc="best")
    style_texts(leg.get_texts(), genotypes)
    sns.despine(ax=ax, offset=10, trim=True)
    fig.savefig(out("plot_morisita_20um.png"), dpi=300)
    plt.close(fig)


# ==================== main analysis ====================

def run_analysis(workbook_path):
    _load_analysis_libs()
    sns.set_theme(style="ticks", context="talk")
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]

    output_dir = os.path.join(os.path.dirname(os.path.abspath(workbook_path)), ANALYSIS_DIR_NAME)
    os.makedirs(output_dir, exist_ok=True)
    out = lambda n: os.path.join(output_dir, n)

    wb = load_workbook(workbook_path, data_only=True)
    genotypes = [s for s in wb.sheetnames if s not in UTILITY_SHEETS]
    print(f"[Info] Genotype sheets: {genotypes}")

    res, fano_rows = {}, []
    for g in genotypes:
        ids, dens, lengths, peaks = read_sheet(wb, g)
        imds, first, nf, fano, mor = compute_spatial_metrics(peaks, COMP_RANGE_20)
        res[g] = dict(ids=ids, densities=dens, lengths=lengths, imds=imds, first_mit=first,
                      norm_freq=nf, fano=fano, morisita=mor)
        print(f"[Info] {g}: {len(ids)} animals, {sum(len(p) for p in peaks)} peaks")
        fano_rows.append({"Genotype": g, "Fano Factor (Mean)": np.nanmean(fano), "Fano Factor (SEM)": sem(fano)})

    fano_df = pd.DataFrame(fano_rows)
    print("\n================== FANO FACTOR ==================")
    print(fano_df.to_string(index=False))

    # ---- summary workbook ----
    with pd.ExcelWriter(out("Genotypes_Summary.xlsx"), engine="openpyxl") as writer:
        rows = []
        for g in genotypes:
            r = res[g]
            rows.append({"Genotype": g, "N animals": len(r["ids"]),
                         "Density Mean": np.nanmean(r["densities"]), "Density SEM": sem(r["densities"]),
                         "Length Mean": np.nanmean(r["lengths"]), "Length SEM": sem(r["lengths"]),
                         "First Mit Dist Mean": np.nanmean(r["first_mit"]), "First Mit Dist SEM": sem(r["first_mit"])})
            width = max((len(x) for x in r["imds"]), default=0)
            padded = [list(x) + [np.nan] * (width - len(x)) for x in r["imds"]]
            pd.DataFrame(np.array(padded).T, columns=r["ids"]).to_excel(
                writer, sheet_name=f"{g}_Step_Distances"[:31], index=False)
        pd.DataFrame(rows).to_excel(writer, sheet_name="Global_Morphology_Summary", index=False)
        def msd(x):
            x = np.asarray(x, dtype=float)
            x = x[~np.isnan(x)]
            return f"{np.mean(x):.2f} \u00b1 {np.std(x, ddof=1):.2f}"
        sd_rows = [{"Genotype": g,
                    "N animals": len(res[g]["ids"]),
                    "Density (per 100 um)": msd(res[g]["densities"]),
                    "Process length (um)": msd(res[g]["lengths"]),
                    "Distance of first mit (um)": msd(res[g]["first_mit"]),
                    "Fano factor": msd(res[g]["fano"])} for g in genotypes]
        pd.DataFrame(sd_rows).to_excel(writer, sheet_name="Mean_SD_Summary", index=False)
        fano_df.to_excel(writer, sheet_name="Fano_Factor", index=False)

    # ---- statistics ----
    print("\n================== STATISTICS ==================")
    metrics = {
        "Mitochondrial Density": "densities",
        "Process Length": "lengths",
        "Distance to First Mit": "first_mit",
    }
    stats_rows, posthoc = [], {}
    for name, key in metrics.items():
        groups = [np.asarray(res[g][key], dtype=float) for g in genotypes]
        groups = [x[~np.isnan(x)] for x in groups]
        print(f"\nMetric: {name}")
        all_normal = True
        for g, x in zip(genotypes, groups):
            if len(x) < 3:
                print(f"  Shapiro-Wilk ({g}): N={len(x)} too small, treating as non-normal")
                all_normal = False
                continue
            p = shapiro(x)[1]
            print(f"  Shapiro-Wilk ({g}): p = {p:.4f} " + ("(normal)" if p > 0.05 else "(NON-NORMAL)"))
            all_normal &= p > 0.05

        if all_normal:
            test, (stat, p) = "One-Way ANOVA", f_oneway(*groups)
            print(f"  One-way ANOVA: F = {stat:.4f}, p = {p:.4f}")
            if p < 0.05:
                print("  Significant -> Tukey HSD p-values:")
                t = tukey_hsd(*groups)
                mat = pd.DataFrame(t.pvalue, index=genotypes, columns=genotypes)
                print(mat.to_string())
                posthoc[(name, "Tukey HSD")] = mat
        else:
            test, (stat, p) = "Kruskal-Wallis", kruskal(*groups)
            print(f"  Kruskal-Wallis: H = {stat:.4f}, p = {p:.4f}")
            if p < 0.05:
                print("  Significant -> Dunn's post-hoc (Holm):")
                flat = pd.DataFrame({"Value": np.concatenate(groups),
                                     "Group": np.repeat(genotypes, [len(x) for x in groups])})
                mat = sp.posthoc_dunn(flat, val_col="Value", group_col="Group", p_adjust="holm")
                mat = mat.loc[genotypes, genotypes]
                print(mat.to_string())
                posthoc[(name, "Dunn")] = mat
            else:
                print("  No significant difference between genotypes.")
        stats_rows.append({"Metric": name, "Test": test, "Statistic": round(stat, 4),
                           "p-value": round(p, 4), "Significant (p<0.05)": "Yes" if p < 0.05 else "No"})

    stats_df = pd.DataFrame(stats_rows)
    save_table_figure(fano_df, "Fano Factor Summary", out("table_fano_factor.png"), italic_first_col=True)
    save_table_figure(stats_df, "Statistical Test Summary", out("table_stats_summary.png"))
    for (name, kind), mat in posthoc.items():
        save_table_figure(mat.reset_index().rename(columns={"index": "Group"}),
                          f"{kind} Post-Hoc ({name})",
                          out(f"table_{kind.split()[0].lower()}_{name.replace(' ', '_').lower()}.png"),
                          italic_first_col=True, italic_header=True)

    pvals = {metrics[name]: mat for (name, _kind), mat in posthoc.items()}
    generate_plots(res, genotypes, output_dir, pvals)
    print(f"\n[Done] Outputs saved in: {output_dir}")


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
        if os.path.isdir(os.path.join(ROOT, d)) and d != ANALYSIS_DIR_NAME
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
