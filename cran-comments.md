## Resubmission

This is a resubmission of 0.3.0. The previous upload was auto-rejected
because two \eqn expressions in simulate_metafrontier.Rd contained
corrupted control characters (tab and backspace), which broke the PDF
manual build. The documentation source has been repaired, and R CMD
check with the PDF manual now passes cleanly. A stray Rplots.pdf test
artifact has also been removed from tests/.

## R CMD check results

0 errors | 0 warnings | 0 notes

(CRAN incoming checks may show 1 NOTE for possibly misspelled words:
QP, metatechnology, poolability, pre (from "pre-fitted"). These are
intentional domain-specific terms.)

## Release summary

metafrontier 0.3.0 (previous CRAN version: 0.2.2).

- New default efficiency estimator: the Battese-Coelli (1988)
  conditional expectation (`estimator = "bc88"`); the JLMS estimator
  remains available, both are stored on every fit, and the change is
  documented as a breaking change in NEWS.md.
- Both identification criteria for the deterministic metafrontier
  (`objective = c("lp", "qp")`; quadprog added to Suggests).
- Backend delegation via `engine = c("internal", "sfaR", "frontier",
  "Benchmarking")`.
- `malmquist_meta()` gains id-based firm matching with strict checks
  and counted infeasibility warnings.
- New `check_convergence()` diagnostic; convergence reporting in
  `print()` and `summary()`.
- DEA additions: FDH, hyperbolic efficiency, user-supplied direction
  vectors, two-stage slack analysis; permutation poolability test for
  DEA fits.
- Bug fixes, including a row-alignment fix for panel SFA
  efficiencies, BC92 decay anchoring on unbalanced panels, and an
  erroneous lower bound in cross-period DEA programs. See NEWS.md.

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
