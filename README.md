# Mitochondrial Distribution Pipeline

A single-file Python tool for measuring how mitochondria are distributed along neuronal processes in fluorescence images. It extracts an intensity profile along a hand-traced ROI, detects candidate mitochondria automatically, lets you correct them in an interactive review panel, stores everything in an Excel workbook, and then runs the statistics and plots.

I built this for my own imaging analysis, so it is a working research tool rather than a polished package.

> **Status:** work in progress. The code is functional and I use it, but it is being refactored.

## What it does

**Features**

- Automated mitochondrial detection
- Interactive review and correction
- Quantitative spatial analysis
- Excel export for correction tracking
- Graphical and tabular analysis output

**Visual Overview**

<table>
  <tr>
    <td align="center" width="60%">
      <img
        src="https://github.com/user-attachments/assets/66f13c79-8009-49c1-964e-39cdbc08d832"
        width="70%"
        alt="MitoSpacing interactive review panel"
      />
      <br>
      <em>Interactive review panel</em>
    </td>
    <td align="center" width="40%">
      <img
        src="https://github.com/user-attachments/assets/fc99c430-8bc5-4cbe-beff-9769a08f8657"
        width="100%"
        alt="MitoSpacing quantitative analysis output"
      />
      <br>
      <em>Quantitative analysis output</em>
    </td>
  </tr>
</table>

**1. Review (interactive)**
- Reads each TIFF image plus its ImageJ/Fiji `.roi` trace and samples the intensity profile along the trace.
- Auto-detects peaks (candidate mitochondria) with `scipy.signal.find_peaks`.
- Opens a two-panel review window: the image with the ROI and arrows marking peaks, and the intensity profile with detected peaks.
  - Click near a green marker to **remove** a peak; click empty space to **add** one.
  - Hovering over the profile highlights the matching position on the image.
  - Add a free-text comment per image.
  - **Save & Close** keeps the image open for later review; **Mark Done & Close** marks it reviewed so future runs skip it.
- Progress is never lost: closing the window saves, and if the workbook is locked (e.g. open in Excel) edits go to a timestamped backup file.
- Every add/remove is recorded in a `corrections_log` sheet, and files that could not be paired are listed in `unmatched_files`.

**2. Analysis**
Once every image is reviewed (or with `--analyze`), the script computes:
- Mitochondrial density (per 100 µm), process length, and distance of the first mitochondrion from the cell body
- Intermitochondrial distance (IMD) distributions, binned in 3 µm steps
- Fano factor per process
- Morisita index across a range of compartment sizes (10 µm and 20 µm steps)
- Group comparisons: Shapiro-Wilk normality check, then one-way ANOVA or Kruskal-Wallis, with Dunn's post-hoc test (Holm-corrected) when applicable

Outputs (saved to `analysis_output/` next to the workbook): summary Excel file, violin/box plots, histograms per genotype, Morisita plots, and table images of the statistics.

## Quantitative analyses

- Mitochondrial density
- Neurite length
- Distance of first mitochondrion to the cell body
- Normalised frequency histogram of inter-mitochondrial distances
- Fano factor 
- Morishita's index of dispersion plotted against bin size

## Expected folder layout

```
ROOT/                     (one date folder)
├── genotype_1/
├── genotype_2/
└── ...
```

Inside each genotype folder, either layout works and is detected automatically:

- **Flat:** image and `.roi` files sit together and are matched by filename (words like `mosaic` / `roi` are ignored when matching).
- **Nested:** one subfolder per animal, each containing one image and one `.roi` file.

Unmatched files are flagged in the workbook rather than silently skipped.

## Installation

Requires Python 3.9+ (tkinter ships with most Python installers).

```bash
pip install -r requirements.txt
```

## Usage

1. Open `mitospacing.py` and edit the config block at the top:
   - `ROOT`: path to your date folder
   - `PX_PER_UM`: your image calibration (pixels per µm)
   - `PEAK_DISTANCE`, `PEAK_PROMINENCE`: peak-detection sensitivity
2. Run the review loop:
   ```bash
   python mitospacing.py
   ```
   When every image is marked done, a dialog offers to run the analysis.
3. To re-run analysis on the current workbook without reviewing anything:
   ```bash
   python mitospacing.py --analyze
   ```

To re-open an image you already marked done, clear its `reviewed` cell in the workbook, save and close the workbook, then rerun the script.

<!--
## Visual Overview

<img width="800" height="450" alt="mitospacing_demo" src="https://github.com/user-attachments/assets/66f13c79-8009-49c1-964e-39cdbc08d832" />
<img width="956" height="474" alt="Screenshot 2026-09-21 133547" src="https://github.com/user-attachments/assets/cade5343-1f08-433f-a15d-2e5d5f2a6f18" />
<img width="582" height="329" alt="Screenshot 2026-09-28 180139" src="https://github.com/user-attachments/assets/fc99c430-8bc5-4cbe-beff-9769a08f8657" />
-->


## Known limitations

- Configuration is edited directly in the script (no command-line arguments or config file yet).
- Bin ranges and compartment sizes for the spatial statistics are hard-coded and tuned to my data.
- Written and tested on Windows.

## Data

No data is included in this repository. Example images and outputs shown above are for illustration only.

## License

MIT License. See [LICENSE](LICENSE) for details.
