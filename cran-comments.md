## R CMD check results

0 errors | 0 warnings | 0 notes

(CRAN incoming checks may show 1 NOTE for new submission and
possibly misspelled words: Battese, Malmquist, Metafrontier, SFA.
These are correct domain-specific terms.)

## Resubmission

Previous submission returned 2 WARNINGs due to missing pre-built
vignettes in inst/doc. This has been fixed — all three vignettes
are now pre-built and included.

## Test environments

* Local: Windows 11, R 4.5.1
* GitHub Actions:
  - Ubuntu 24.04, R-devel
  - Ubuntu 24.04, R 4.5.1 (release)
  - Ubuntu 24.04, R 4.4.x (oldrel-1)
  - Windows Server 2022, R 4.5.1 (release)
  - macOS 14, R 4.5.1 (release)

## Notes

The CRAN incoming checks report two DOIs as "404 Not Found":

* `10.1023/B:PROD.0000012454.06094.29` (Battese, Rao, and O'Donnell, 2004)
* `10.1007/s11123-014-0402-2` (Huang, Huang, and Liu, 2014)

Both DOIs are valid and resolve correctly in a web browser. The 404
responses are caused by the doi.org resolver returning errors on
programmatic HEAD requests for certain older Springer/Kluwer DOIs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
