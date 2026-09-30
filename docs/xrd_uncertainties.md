# XRD intensity columns and fit statistics

Version 1.2.3 preserves input intensity uncertainties without instrument-preset
rescaling. For a named Synergy/CrysAlisPro CSV, the fitter selects columns by
their headers rather than their positions.

| Header | Meaning in the powder export | Use in the fit |
| --- | --- | --- |
| `2thetadeg` | Diffraction angle 2θ, in degrees | Horizontal coordinate |
| `intx` | Processed intensity for an angular bin | Observed intensity |
| `d-value` | Bragg-equivalent spacing in Å, d = λ / (2 sin θ) | Not used in a fit against 2θ |
| `sigx` | Exported intensity uncertainty, interpreted as one standard deviation | Inverse-variance weight: w = 1 / sigx² |
| `count` | A contribution/coverage counter in this export; exact vendor aggregation is not independently documented here | Exclude bins with count ≤ 0; positive count is not an extra weight |

CrysAlisPro extracts one-dimensional powder profiles from detector images and
can apply different integration modes and corrections. Therefore, `intx` should
not automatically be described as raw photon counts or counts per second.
The `count` values behave as coverage/contributions, but the precise distinction
between pixels and accumulated contributions across images requires the vendor's
export specification. It is not interchangeable with `intx`.

Use the supplied `sigx` in the same intensity units as `intx`. Do not divide it
again by sqrt(count), replace it with sqrt(intx), or multiply it by an arbitrary
factor to obtain a preferred goodness of fit. If no uncertainty column is
available, the existing sqrt(intensity) fallback is only a counting-statistics
estimate; normalized intensity generally needs additional acquisition information.
Measured rows with invalid supplied uncertainties are rejected by the named CSV
reader. Intermediate GSAS XYE files retain floating-point precision, including
small positive uncertainties.

The fit minimizes a weighted residual sum, Σ[(Iobs − Icalc) / sigx]². In the
usual unrestrained case, reduced chi-squared divides that sum by N − p, where
N is the number of included observations and p the number of refined independent
parameters. GOF is its square root. The toolkit reports the native GSAS-II
statistics, which also handle the refinement's own exclusions and restraints.

GOF near one is expected only when the model and statistical assumptions are
appropriate. A larger value can reflect model mismatch, incomplete uncertainty
estimates or correlations in processed data. A small value is not proof of a
better structural model. Inspect residuals, physical parameters and processing
assumptions rather than adjusting errors to force agreement.

## Earlier versions

The old Synergy-S preset specified a fivefold sigma adjustment. That was an
empirical choice for another dataset, not a general detector calibration.
It was also written to the wrong GSAS-II histogram-tree location, leaving native
weights unchanged while the displayed GOF was reduced. A separate positional
CSV import error could use `d-value` as sigma, changing relative weights and
potentially the fitted parameters. Version 1.2.3 removes the preset adjustment,
uses native statistics and fixes named CSV import. Refit affected scans; changing
only a displayed GOF cannot repair a fit performed with the wrong weights.

## References and limits of the definitions

- [Rigaku application note SMX027, page 6](https://rigaku.com/hubfs/2024%20Rigaku%20Global%20Site/Resource%20Hub/Applications%20Library/PDF/145539922978.pdf?hsLang=en) describes 2D-to-1D powder extraction and processing choices. It does not define these exact CSV headers.
- [CrysAlisPro user manual](https://www.agilent.com/cs/library/usermanuals/public/CrysAlis_Pro_User_Manual.pdf), section 6.2, describes powder extraction and XYE/GSAS export.
- [Toby (2024), Journal of Applied Crystallography](https://doi.org/10.1107/S1600576723011032), section 3, describes intensity uncertainties, weighting and GOF.
