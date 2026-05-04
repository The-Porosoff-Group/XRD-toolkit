# XRD Fitting Toolkit

Standalone XRD front end for the Catalysis Data Toolkit.

This launcher opens directly into the XRD fitting workflow and hides the unfinished GC/TGA/BET navigation. It uses the same backend as the full toolkit, including Materials Project CIF import, CIF caching, GSAS-II refinement, saved presets, phase-card controls, baseline comparison, and fit warnings.

## Quick Start

1. Add your Materials Project API key:

   ```bat
   copy ..\config.yaml.example ..\config.yaml
   notepad ..\config.yaml
   ```

2. Launch the standalone XRD GUI:

   ```bat
   run_xrd_toolkit.bat
   ```

3. The browser opens at:

   ```text
   http://localhost:5000/xrd
   ```

The first run may take several minutes while Python dependencies and GSAS-II are installed.

## What This Includes

- XRD data upload and live preview
- Materials Project phase search by elements, formula, or keyword
- Manual CIF upload
- CIF validation and cached CIF handling
- Preview tick generation from imported phases
- GSAS-II refinement
- Instrument parameter file support (`.instprm`)
- Built-in WC/W2C Synergy-S production preset
- Saved user presets in `xrd_refinement_presets.json`
- Per-phase refinement cards for size, microstrain, March-Dollase preferred orientation, and W2C-style uniform-cell diagnostics
- Fit quality statistics, phase fractions, uncertainty notes, FWHM reference peaks, crystallite sizes, preferred-orientation results, and cell-change percentages
- Baseline comparison and warnings for overfitting-prone settings

## Recommended Fitting Workflow

1. Load an XRD pattern.
2. Set the wavelength and 2-theta range.
3. Search Materials Project and queue the expected phases.
4. Run a constrained baseline fit first.
5. Save or mark that baseline in the GUI.
6. Add one model freedom at a time, such as cell refinement, preferred orientation, Uiso, size, or microstrain.
7. Compare the new fit against the baseline.
8. Keep extra parameters only when they improve the expected residual features and the fitted values remain physically reasonable.
9. Save the validated recipe as a preset for related samples.

For the WC/W2C workflow, the shipped preset fixes WC [001] March-Dollase preferred orientation near `0.905` because that value was established from repeated comparison fits. It is a production prior for that specific recipe, not a universal WC constant.

## Important Files

From the toolkit root:

```text
app.py                         Flask backend and routes
templates/xrd_toolkit/index.html
                               XRD-only front end
xrd_toolkit/run_xrd_toolkit.bat
                               Standalone Windows launcher
xrd_toolkit/README.md          This standalone module README
modules/xrd/                   XRD parsing, CIF, crystallography, and GSAS-II code
fixtures/                      Canonical CIF fixtures
config.yaml                    Local API-key config, not committed
xrd_refinement_presets.json    Local saved presets, not required for a clean install
```

## Full Toolkit Compatibility

The original full toolkit still runs from:

```bat
..\run.bat
```

and opens:

```text
http://localhost:5000
```

This XRD launcher only changes the startup URL to `/xrd`. It does not remove or disable the full application.
