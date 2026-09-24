# Release notes

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
