# Release notes

## 1.2.2 — 2026-09-24

- Make the XRD workbook's **Fit Parameters** tab a readable recipe for repeating
  the fit in the GUI: original phase-card and CIF identities, instrument and
  file identities, scan range, background settings, every checkbox and entered
  value, and software versions.
- Keep detailed native GSAS-II parameters, arrays, and the fitted model in the
  companion `.gpx` project instead of dumping them as JSON into the worksheet.

## 1.2.1 — 2026-09-24

- Restore the instrument selector beside Run Refinement with three shared
  instruments and **None / calibration**.
- Start Si 640g calibration from fresh parameters by selecting None; show only
  its geometry control, with the Cu doublet chosen automatically for Cu data.
- Remove local profile naming, radiation-spectrum selection, and the separate
  calibration checkbox from setup. Calibration results retain file downloads.
- Keep optional `.instprm` overrides in collapsed GSAS controls for sample fits;
  existing local profile files and APIs remain available.
- Move legend placement into the results' Plot settings and simplify Scan
  Settings in both the main toolkit and standalone XRD interface.

- Show elapsed fitting time instead of cycling through simulated progress stages.

## 1.2.0 — 2026-09-24

- Bundle the benchtop Cu flat-plate Si 640g calibration and show it in the
  instrument dropdown in both the full toolkit and standalone XRD interface.
- Add generic geometry presets, .instprm upload, and named local instruments
  that persist across restarts without overwriting bundled calibrations.
- Calibrate Si 640g using its built-in certified cell, explicit geometry and
  radiation spectrum, fixed sample broadening, and independent peak intensities.
- Apply GSAS-II cycle limits correctly, compare width terms again after
  asymmetry is fitted, and support a validated joint X/Y broadening model.
- Initialize all geometry-specific sample parameters for single wavelengths
  and doublets; preserve the correct Gaussian-only checkpoint ranking.
- Provide candidate downloads, calibration reports, and explicit save controls.
- Display the same release version in both interfaces.
