# sobolnp 0.2.0

## Bug fixes

* Request standard errors explicitly in `npreg()` (`se = TRUE`). np 0.80-1
  made them opt-in, which broke the residual bootstrap.
* The residual bootstrap now resamples the residuals of each variable
  separately. Previously all variables were bootstrapped with the residuals
  of the first variable, and `sobolnp()` failed with two input variables.
  **Bootstrap estimates (`Sboot`, `bwsboot`) differ from version 0.1.0.**
* A `bandwidth` vector of length `ncol(X)` is now used one value per
  variable instead of failing.
* The package no longer masks `base::print()` and `graphics::plot()`.

## Changes

* **Improved bootstrap bandwidth selection.** The bootstrap responses are built
  from residuals standardized with a kernel estimate of their scale. For each
  candidate bandwidth, the bootstrap curves give a bias-corrected curve
  (twice the CV curve minus the mean bootstrap curve). The bandwidth that
  keeps the bootstrap indices closest to the CV index is searched between
  `bws / 1000` and `bws`. In simulations with the default second-order
  kernel, `Sboot` reduces the downward bias of the CV estimate. Results
  differ from version 0.1.0, and the minqa dependency is no longer needed.
* `sobolnp()` always returns a `sobolnp` object, also with
  `bootstrap = FALSE` (then `Sboot` and `bwsboot` are `NULL`).
* `sobolnp()` validates its inputs, accepts data frames and names unnamed
  variables `X1, ..., Xp`. Negative bandwidths raise an error instead of
  returning `NA`.
* New argument `verbose` to silence progress messages and bars.
* `plot()` draws a grouped bar plot, accepts `ylim`, `xlab`, `ylab`, `col`
  and `legend.pos`, and returns its input invisibly.
* `print()` gains a `digits` argument and returns its input invisibly.

# sobolnp 0.1.0

Initial version based in the paper

Solís, Maikol. "Nonparametric estimation of the first order Sobol indices with bootstrap bandwidth." \emph{arXiv preprint arXiv:1803.03333} (2018).
