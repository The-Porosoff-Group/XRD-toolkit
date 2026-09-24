# XRD Fitting Toolkit

**Release 1.2.0** — shared benchtop calibration and local instrument profiles.

Standalone XRD fitting interface from the Catalysis Data Toolkit.

This version opens directly into the XRD workflow and hides unfinished modules. It is intended as a local desktop web app for importing CIFs, previewing XRD tick patterns, and running GSAS-II refinements from a browser-based GUI.

## Quick Start

1. Download or clone this repository.
2. Add your Materials Project API key:

   ```bat
   copy config.yaml.example config.yaml
   notepad config.yaml
   ```

3. Start the XRD-only GUI:

   ```bat
   xrd_toolkit\run_xrd_toolkit.bat
   ```

4. Open, if the browser does not open automatically:

   ```text
   http://localhost:5000/xrd
   ```

The first launch creates a local Python environment and installs dependencies. GSAS-II installation can take several minutes.

## Instrument calibration and local profiles

**Instrument Settings → Instrument / calibration profile** now offers generic
flat-plate (Bragg–Brentano) and capillary/transmission starting profiles, named
bundled profiles, and **Upload .instprm…**. Generic profiles do not load another
instrument's measured calibration. Unidentified scans default to the generic
flat-plate geometry; verify the geometry before fitting.

The **Benchtop Cu — flat plate (Si 640g)** entry loads the bundled
`benchtop_Cu_Si640g.instprm` automatically. It is included in Git for other users
of that same instrument. The selected 20–90° Cu-doublet calibration achieved
Rwp 8.502%; see [settings, comparisons, and limits](docs/benchtop_calibration.md).
The instrument make/model is not recorded, so this is not a generic calibration
for other benchtop machines.

To calibrate from **NIST Si 640g**:

1. Upload the standard measured with the same optics and configuration as your
   samples. Select its geometry and an angular range with several Si peaks.
2. Set the radiation spectrum: Cu Kα1/Kα2 for a normal Cu doublet, or **Single
   wavelength / Kα1 only** for monochromated or Kα2-stripped data. Check beam
   polarization against the optics. Auto spectrum uses a Cu doublet near 1.54 Å
   and a single wavelength for other sources.
3. Tick **Calibrate instrument using NIST Si 640g** and run **GSAS-II Refinement**.
   The certified Si cell (`a = 5.431109 Å`) is built in; no database/API key or
   selected sample phase is required. Selected sample phases are ignored.
4. Inspect the fit and warnings, then download the candidate `.instprm`, report,
   and GSAS-II project. Parameter checks cover the fitted angular range; they
   are not a fit-quality certification. Chi-square requires valid uncertainties,
   which cannot be inferred from counts-per-second data without counting times.
5. Enter a descriptive name and click **Save as local instrument**. The saved
   profile appears immediately in the dropdown and calibration mode turns off.
   Test it with a repeat standard scan before relying on sample broadening.

The calibrator explicitly applies the selected geometry. It starts from fresh
profile values, fixes the certified cell and sample broadening, initializes the
reflection list, then extracts Si peak intensities independently with Le Bail
while fitting instrument widths and Zero in stages. This avoids forcing texture
or slit-dependent intensity ratios into width parameters. Cu doublet intensity
ratio is fixed at 0.5. A mounting-dependent position error can still affect Zero;
check peak positions and repeatability. Instrument/optics-specific asymmetry and
absorption can require further work in GSAS-II.

To reuse an existing calibration, select **Upload .instprm…**, choose its
geometry, select a single-bank GSAS-II `Type:PXC` file, and either fit directly
or enter a name and click **Save uploaded file locally**. The file supplies the
wavelength and spectrum during sample fitting. Legacy `.prm`/`.inst`, TOF, and
multi-bank files must first be converted/exported as a supported `.instprm`.
Other calibration standards can be fitted externally and imported this way.

Local profiles and geometry metadata are stored in the git-ignored
`local_instruments/` folder. They survive restarts, remain on that computer, and
are not bundled in software updates. Each save creates a new profile; existing
files are not overwritten. Calibration candidates and reports are kept together
in that run's results folder. A failed candidate can be downloaded for diagnosis
but cannot be saved from the calibration result as a local instrument.

## Main Features

- XRD file upload with live preview (`.dat`, `.xy`, `.xye`, `.csv`, `.txt`, `.xlsx`)
- Materials Project phase search
- Manual CIF upload
- Versioned CIF caching and pre-refinement validation
- Correct preview tick generation from imported phases
- Conventional-setting Materials Project import with cell/coordinates kept together
- GSAS-II refinement backend
- `.instprm` instrument parameter support
- Built-in WC/W2C Synergy-S production preset
- Saved user presets
- GSAS HAP, Scherrer-equivalent, or combined size reporting with an explicit K
- Light/dark figure exports with an optional title; single-phase legends omit
  the normalized weight percentage
- Legend position selector, including automatic and outside-right placement;
  finished figures can be updated without refitting
- **Fit Parameters** workbook tab with input scans/uncertainties, fit settings,
  initial/final native GSAS-II parameters and flags, refinement stages,
  CIF/instrument inputs, and software versions; companion `.gpx` project download
- Per-phase controls for:
  - unit-cell refinement
  - crystallite size
  - microstrain
  - March-Dollase preferred orientation
- Fit convergence/correlation warnings and baseline comparison
- Outputs for phase fraction, uncertainty notes, exact GSAS profile FWHM,
  crystallite size, microstrain, preferred-orientation value, and cell-change
  percentages relative to the exact prepared GSAS CIF

**Wt%** reports GSAS-II mass fractions normalized over the modeled crystalline
phases. It does not include amorphous material or unmodeled phases. **Diffraction
area (%)** is a separate intensity diagnostic, not a weight percentage. Le Bail
and legacy in-house Rietveld fits report diffraction area only; their Wt% is
unavailable. Refit and regenerate older results to replace percentages that were
previously mislabeled as Wt%.

Figure exports use Arial for text, subscripts, and crystallographic overbars.
Install Arial on the computer generating the figures; the font is not bundled.
If unavailable, exports warn and use Liberation Sans, then DejaVu Sans. Existing
PNGs keep their original typography until regenerated.

The **Fit Parameters** tab stores values as JSON text to retain their precision.
Join numbered parts before decoding a long value. Open the companion `.gpx`
project in GSAS-II to inspect the fitted model; use matching software versions
to reproduce a fit. Older workbooks require a new fit to capture the native state.
Batch runs also accept `--legend-location "outside right"`.

## Recommended Workflow

Raw scans, generated results, and `xrd_refinement_presets.json` are local user
files and are excluded from Git. Keep measurements in `data/` or `uploads/`
and fit outputs in `results/`. Canonical `fixtures/` and the bundled
instrument reference profiles remain versioned. Ignore rules do not remove
files from older Git commits.

1. Upload the measured XRD pattern.
2. Set wavelength and 2-theta range.
3. Search Materials Project or upload CIFs for the expected phases.
4. Confirm the phase card reports `CIF ok`, the intended space group,
   conventional cell, asymmetric-site count, and plausible preview ticks.
5. Run a constrained baseline fit with phase Cell, Size, Mustrain, PO, Uiso,
   and atom positions fixed.
6. Mark/save the baseline in the GUI.
7. Add one refinement freedom at a time. Refine cell before sample
   broadening, and use a measured `.instprm` whenever size or strain matters.
8. Keep a freedom only when the refinement converges, the residual improves
   in the expected peaks, and the parameter remains physically plausible.
9. Save a validated preset for related samples.

Do not keep extra fit freedoms just because Rwp improves. Preferred orientation, Uiso, size, microstrain, and atom-position refinement can all improve the statistic while also changing phase fractions or absorbing model error.

The toolkit rejects a Materials Project phase before refinement if GSAS-II
imports a different space group, cell setting, or asymmetric-site model than
the prepared CIF. Raw P1/full-cell fixtures are retained only as regression
inputs and cannot silently replace normal-import CIFs. Cached MP structures use
versioned normal-import keys so older transformed CIFs cannot poison a run.

## Mo2C M1 Acceptance Recipe

For the included `mp-1552` orthorhombic Mo2C model and a calibrated SmartLab
instrument profile:

1. Fit `20-60` degrees 2-theta with Cu K-alpha (`1.54056` A).
2. Add `mp-1552` and verify Pbcn (No. 60), two asymmetric sites, and cell near
   `a=4.7285`, `b=6.0526`, `c=5.2098` A.
3. Run the constrained baseline with all phase-card freedoms off.
4. Enable Cell for Mo2C and the global `Refine cell too` gate; rerun.
5. Enable Size for Mo2C and rerun with Mustrain still off.

The current acceptance data reaches about `Rwp=5.93%`, a phase-profile FWHM of
`0.445 degrees` for the (121) reference near `39.58 degrees`, and a GSAS HAP
size near `22.4 nm`. Enabling Mustrain changes Rwp by only about `0.01%` and
produces approximately `99%` size/strain correlation, so strain should be
rejected for this example.

## Terminal Batch Workflow

The same refinement backend can run without the GUI.

Use the GUI to build and save a validated recipe, then run many samples from the terminal:

```bat
xrd_toolkit\run_xrd_batch.bat ^
  --patterns data\*.xy ^
  --preset "WC/W2C Synergy-S production" ^
  --cif-dir cifs\wc_w2c ^
  --out results\wc_w2c_batch
```

The batch wrapper performs the same first-run environment setup as `run.bat`, then calls `scripts\xrd_batch.py`.

You can provide phases in any of these ways:

- GUI-saved preset phases from `xrd_refinement_presets.json`
- local CIF folder:

  ```bat
  --cif-dir cifs\wc_w2c
  ```

- one or more explicit CIFs:

  ```bat
  --cif cifs\WC.cif --cif cifs\W2C.cif
  ```

- Materials Project ids:

  ```bat
  --mp-ids mp-2034 mp-your-wc-id
  ```

- a standalone recipe JSON:

  ```bat
  --recipe recipes\wc_w2c_synergy_s_batch.json
  ```

Batch outputs include one folder per sample plus aggregate summaries:

```text
results/wc_w2c_batch/
  resolved_batch_recipe.json
  batch_summary.json
  batch_phase_summary.csv
  sample_1/
    xrd_refinement.png
    xrd_summary.xlsx
    summary.json
```

To pre-generate CIFs from Materials Project:

```bat
xrd_toolkit\fetch_cifs.bat --mp-ids mp-2034 mp-your-wc-id --out-dir cifs\wc_w2c
```

## WC/W2C Preset

The built-in WC/W2C Synergy-S preset uses a fixed WC [001] March-Dollase preferred-orientation value near `0.905`. That value came from a comparison workflow and is meant as a production prior for this specific recipe. It is not a universal WC constant.

## Important Files

```text
app.py                         Flask backend and routes
run.bat                        Full toolkit launcher
xrd_toolkit/run_xrd_toolkit.bat
                               XRD-only launcher
xrd_toolkit/run_xrd_batch.bat  XRD batch launcher
xrd_toolkit/fetch_cifs.bat     Materials Project CIF fetch launcher
templates/xrd_toolkit/index.html
                               XRD-only GUI
scripts/xrd_batch.py           Terminal batch refinement wrapper
scripts/fetch_cifs.py          Materials Project CIF fetch helper
recipes/                       Terminal recipe JSON files
modules/xrd/                   XRD, CIF, crystallography, and GSAS-II code
fixtures/                      Canonical CIF fixtures
config.yaml.example            API-key template
```

## Notes on GitHub Pages

This toolkit cannot run directly as a static GitHub Pages site because it depends on Python, Flask, GSAS-II, local file uploads, and local refinement outputs. GitHub can host the source code and documentation, but users run the app locally with the launcher.
