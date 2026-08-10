## Release summary

metafrontier 0.3.1 (previous CRAN version: 0.3.0).

This is a patch release addressing the failure reported by the CRAN
additional check without long doubles (noLD): the test
"bc88 is computed in all distribution branches" failed for the
heteroscedastic exponential model.

Cause: the BC88/BC92 conditional efficiencies were evaluated as a
natural-scale ratio of two `pnorm()` values. In degenerate fits with
`sigma_u` near zero, both values are subnormal, and without
extended-precision long doubles the ratio can return 0, `Inf`, or
values above one.

Fix: the efficiencies are now evaluated on the log scale
(`pnorm(log.p = TRUE)`), with an asymptotic closed form in the far
left tail where the log-scale difference itself loses precision. On
standard platforms results are unchanged to near machine precision.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.5.1
* GitHub Actions:
  - Ubuntu 24.04, R-devel
  - Ubuntu 24.04, R 4.5.1 (release)
  - Ubuntu 24.04, R 4.4.x (oldrel-1)
  - Windows Server 2022, R 4.5.1 (release)
  - macOS 14, R 4.5.1 (release)

## Notes

The CRAN incoming checks may report two DOIs as "404 Not Found":

* `10.1023/B:PROD.0000012454.06094.29` (Battese, Rao, and O'Donnell, 2004)
* `10.1007/s11123-014-0402-2` (Huang, Huang, and Liu, 2014)

Both DOIs are valid and resolve correctly in a web browser. The 404
responses are caused by the doi.org resolver returning errors on
programmatic HEAD requests for certain older Springer/Kluwer DOIs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
