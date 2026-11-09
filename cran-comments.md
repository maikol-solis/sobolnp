# Version 0.2.0

## Submission

This release fixes check failures reported by the np maintainer
(J. S. Racine) against np 0.80-1, where standard-error computation in
`npreg()` became opt-in. `sobolnp` now requests `se = TRUE` explicitly when
the uncertainty vector is needed.

It also fixes a bug in the residual bootstrap, stops masking the base
`print()` and `plot()` generics, and adds input validation and testthat
tests.

## Test environments

* local macOS, R 4.6.2, np 0.80-1

## R CMD check results

0 errors | 0 warnings | 0 notes
