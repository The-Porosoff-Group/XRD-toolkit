# Benchtop Cu flat-plate calibration — v1.2.0

Select **Benchtop Cu — flat plate (Si 640g)** in Instrument Settings in either
the Catalysis Data Toolkit or the standalone XRD tool. The bundled file is
[`benchtop_Cu_Si640g.instprm`](../benchtop_Cu_Si640g.instprm). It is included in
Git, so pulling this release also installs the dropdown entry. No upload is
needed. A separately saved local instrument remains local to its installation.

This empirical calibration belongs to the benchtop instrument that produced
the supplied `Standard.txt` scan. Its make/model and detailed optics were not
recorded. It is **not a universal calibration for all benchtop instruments**.
Use the same instrument, geometry, optics, and data processing as the standard.

## Selected settings

Calibrated on 2026-09-24 using NIST Si 640g, with the certified cubic cell fixed
at **a = 5.431109 Å**. The certified reference temperature is 22.5 °C; the scan's
measurement temperature was not supplied. NIST describes this standard's
sample broadening as negligible for laboratory instruments; the GSAS-II model
uses fixed 10 µm size and zero microstrain.
[NIST certification](https://doi.org/10.6028/NIST.SP.260-245).

| Setting | Selected value |
|---|---|
| Geometry | Flat plate / Bragg-Brentano |
| Fitted range | 20–90° 2θ; six Si reflection families |
| Radiation | Cu Kα1 = 1.540593 Å; Kα2 = 1.544414 Å |
| Kα2/Kα1 ratio | 0.5, fixed |
| Polarization | 0.5, assumed; peak intensities independently extracted |
| Background | Six Chebyshev coefficients |
| Instrument parameters fitted | W, X, Y, Zero, SH/L |
| U, V, Z | Fixed at zero in the selected profile |
| Sample height and transparency | Fixed at zero for export |
| Final Rwp / Rp | 8.502% / 6.231% |

The selected exported values are Zero = −0.04099°, W = 4.6106,
X = 2.1253, Y = 6.8818, and SH/L = 0.06399. U/V/W use GSAS-II's Gaussian
variance convention in centidegrees squared; X/Y use its Lorentzian convention.

## Model comparisons

All comparisons below used the same 20–90° range. Reflection intensities were
extracted independently while the cell and sample broadening remained fixed.
Rwp is a measure of agreement with this scan, not independent validation.

| Comparison | Rwp (%) | Decision |
|---|---:|---|
| Single wavelength, six background terms | 21.38 | Reject: resolved high-angle doublets are not represented |
| Doublet, W + X, fixed small asymmetry | 13.61 | Reject: asymmetry mismatch |
| Doublet, W + X with fitted asymmetry | 10.50 | Baseline |
| Doublet, U + W + X with fitted asymmetry | 8.57 | Valid alternative |
| **Doublet, W + X + Y with fitted asymmetry** | **8.50** | **Selected** |
| Add U to W + X + Y | 8.49 | Improvement too small to justify another width term |
| U + V + W + X | 8.24 | Reject: negative Gaussian variance within the fitted range |
| Float the doublet ratio | 8.20 | Reject for this export: fitted ratio ~0.385 conflicts with high-angle pairs near 0.5 |
| Also refine sample height | ~7.9 | Exclude mounting-dependent correction from the shared instrument file |

Increasing the background to 10 or 14 terms barely improved the W + X
baseline. Six terms were retained. The calibrator now revisits width terms
after fitting asymmetry and compares joint X/Y models with physical-range and
correlation checks. Instrument parameters stay fixed when this measured profile
is loaded for ordinary sample fitting.

## Validation and limits

- The selected parameters are finite, have nonnegative X/Y, and predict positive
  Gaussian variance throughout 20–90°. More complex invalid models were rejected.
- High-angle peak separations and independently fitted component ratios support
  ordinary Cu Kα1/Kα2 data. They do not establish the instrument's make/model.
- Zero can include specimen mounting error from this single scan. Check sample
  positioning and, when appropriate, refine the sample-height correction on
  subsequent samples rather than changing the measured instrument widths.
- Residual peak-shape mismatch remains. This is the best-supported tested
  transferable model, not proof of a unique or perfect instrument function.
- Do not extrapolate beyond the calibrated range without checking another
  standard. A repeat or independent standard scan was not provided.
- The source is in counts per second without counting times. Approximate
  square-root-intensity weights were used; chi-square is not a calibrated
  counting-statistics goodness-of-fit test.

The original measurement and native fitting projects remain local. The shared
profile and this settings/selection record are versioned in both repositories.
