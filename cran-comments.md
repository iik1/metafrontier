## Release summary

metafrontier 0.4.0 (previous CRAN version: 0.3.1).

New features: cost frontiers, sample-selection (Greene, 2010) and
spatial-autoregressive group frontiers, and transitive Fare-Primont and
Lowe productivity indexes (`tfp_index()`). Bug fixes, among them a
numerical failure of exponential SFA fits and an incorrect technology
gap ratio for the stochastic metafrontier. See NEWS.md for details.

New suggested packages: productivity, Matrix and spdep (used in tests
and, for Matrix and spdep, optionally for spatial weights input).

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.5.1
* GitHub Actions:
  - Ubuntu 24.04, R-devel
  - Ubuntu 24.04, R release
  - Ubuntu 24.04, R oldrel-1
  - Windows Server 2022, R release
  - macOS 14, R release

## Notes

The CRAN incoming checks may report two DOIs as "404 Not Found":

* `10.1023/B:PROD.0000012454.06094.29` (Battese, Rao, and O'Donnell, 2004)
* `10.1007/s11123-014-0402-2` (Huang, Huang, and Liu, 2014)

Both DOIs are valid and resolve correctly in a web browser. The 404
responses are caused by the doi.org resolver returning errors on
programmatic HEAD requests for certain older Springer/Kluwer DOIs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
