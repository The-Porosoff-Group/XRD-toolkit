# XRD Fitting Toolkit

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
- Per-phase controls for:
  - unit-cell refinement
  - crystallite size
  - microstrain
  - March-Dollase preferred orientation
- Fit convergence/correlation warnings and baseline comparison
- Outputs for phase fraction, uncertainty notes, exact GSAS profile FWHM,
  crystallite size, microstrain, preferred-orientation value, and cell-change
  percentages relative to the exact prepared GSAS CIF

## Recommended Workflow

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
