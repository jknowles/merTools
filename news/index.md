# Changelog

## merTools 1.0.0

This is the 1.0 long-term-support release. It resolves the remaining
open issues, fixes a correctness bug in
[`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
for nested random effects, repairs and extends the
[`shinyMer()`](https://jknowles.github.io/merTools/reference/shinyMer.md)
explorer, adds a new
[`plotREimpact()`](https://jknowles.github.io/merTools/reference/plotREimpact.md)
visualization, refreshes the documentation, and tidies the test suite
for low-maintenance, long-term use.

### New Features

- **New
  [`plotREimpact()`](https://jknowles.github.io/merTools/reference/plotREimpact.md)
  function for visualizing
  [`REimpact()`](https://jknowles.github.io/merTools/reference/REimpact.md)
  output ([\#84](https://github.com/jknowles/merTools/issues/84),
  [\#85](https://github.com/jknowles/merTools/issues/85)).** Plots the
  weighted-average fitted value for each expected-rank bin of a grouping
  factor, with confidence intervals, faceted by case. Passing a *named
  list* of
  [`REimpact()`](https://jknowles.github.io/merTools/reference/REimpact.md)
  results overlays them on shared axes so the influence of different
  grouping factors (or the same factor across models) can be compared
  directly — previously this required hand-assembling the data frames.
  Uses a clean
  [`theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  look.
- **[`plotFEsim()`](https://jknowles.github.io/merTools/reference/plotFEsim.md)
  now highlights significant fixed effects
  ([\#85](https://github.com/jknowles/merTools/issues/85)).** Terms
  whose interval excludes the null line are drawn in solid black and the
  rest in grey, matching the convention already used by
  [`plotREsim()`](https://jknowles.github.io/merTools/reference/plotREsim.md),
  on a cleaner
  [`theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  canvas with labelled axes.
- **[`shinyMer()`](https://jknowles.github.io/merTools/reference/shinyMer.md)
  gained a “Model Summary” tab
  ([\#78](https://github.com/jknowles/merTools/issues/78)).** The
  interactive explorer now opens a dedicated tab summarizing the fitted
  model: the
  [`modelInfo()`](https://jknowles.github.io/merTools/reference/modelInfo.md)
  overview (observations, grouping factors, AIC, residual sigma), the
  original model call, the fixed-effect estimates, the random effect
  variances/correlations (`VarCorr`), and the number of levels per
  grouping factor (`ngrps`).
- **[`shinyMer()`](https://jknowles.github.io/merTools/reference/shinyMer.md)
  can now restrict the Random / Average draw to a subset of the data
  ([\#32](https://github.com/jknowles/merTools/issues/32)).** When the
  “Random Obs” or “Average Obs” scenario is selected, a sidebar control
  lets you pick a model-frame variable and value to filter on; the
  chosen case is drawn from that subset (via the existing
  `draw(..., varList = )` machinery), enabling much richer on-the-fly
  exploration of specific cases.

### Documentation

- **New “Exploring Contextual Effects with merTools” vignette
  ([\#136](https://github.com/jknowles/merTools/issues/136)).** A fresh,
  worked example using the bundled `hsb` data (student SES
  vs. school-mean SES) that separates within-group and contextual
  effects, mirroring the easystats `modelbased` “effects in context”
  article (cited). It demonstrates
  [`FEsim()`](https://jknowles.github.io/merTools/reference/FEsim.md)/
  [`plotFEsim()`](https://jknowles.github.io/merTools/reference/plotFEsim.md),
  [`wiggle()`](https://jknowles.github.io/merTools/reference/wiggle.md) +
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md),
  and the new
  [`plotREimpact()`](https://jknowles.github.io/merTools/reference/plotREimpact.md)
  on a single, self-contained model.
- **Rebuilt the “Prediction Intervals” vignette and fixed its formatting
  ([\#116](https://github.com/jknowles/merTools/issues/116)).** Several
  section headers were missing the space after `#` (so they rendered as
  body text) and two sub-sections were both numbered “Step 3c”; these
  are corrected. The
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)-vs-[`bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html)
  comparison figure, which no longer matched the surrounding text, has
  been regenerated with the current code and now shows the methods
  agreeing closely. A latent
  [`display()`](https://rdrr.io/pkg/arm/man/display.html) call (from
  `arm`, which was not attached) was qualified to
  [`arm::display()`](https://rdrr.io/pkg/arm/man/display.html).
- **Improved package citation and acknowledged influences
  ([\#137](https://github.com/jknowles/merTools/issues/137)).**
  `inst/CITATION` now derives its version and year from the package
  metadata (rather than hard-coding a stale version), lists all package
  authors, and a citation footer points users to the methodological
  foundations — Gelman and Hill (2007), the `arm` package’s
  [`sim()`](https://rdrr.io/pkg/arm/man/sim.html), and `lme4` (Bates et
  al. 2015). The package-level help page
  ([`?merTools`](https://jknowles.github.io/merTools/reference/merTools-package.md))
  gained matching “Influences and acknowledgements” and “References”
  sections.

### Bug Fixes

- **Fixed two
  [`shinyMer()`](https://jknowles.github.io/merTools/reference/shinyMer.md)
  defects in the “Substantive Effect” tab.** The fixed-effect impact
  (“wiggle”) plot crashed for any non-numeric fixed effect
  (`object 'newvals' not found`), and the numeric/factor branch used
  `class(x) %in% ...` inside `if()`, which errors with “the condition
  has length \> 1” for multi-class objects (e.g. ordered factors) on R
  \>= 4.2. Both are corrected
  ([`is.numeric()`](https://rdrr.io/r/base/numeric.html) dispatch, and
  factor/character values now wiggle across their observed levels), and
  the per-case faceting index is computed correctly. The app’s data
  table was also migrated from the deprecated
  [`shiny::renderDataTable()`](https://rdrr.io/pkg/shiny/man/renderDataTable.html)
  to
  [`DT::renderDT()`](https://rdrr.io/pkg/DT/man/dataTableOutput.html),
  the server now declares a `session` argument, and the
  substantive-effect interval uses a correct two-sided width.

- **Fixed non-reproducible
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  for partially-observed nested / interaction random effects
  ([\#124](https://github.com/jknowles/merTools/issues/124)).** For a
  model such as `y ~ 1 + (1 | a / b)`, a prediction frame mixing
  observed and unobserved levels of the interaction grouping factor
  (`a:b`) returned seed-dependent point estimates, and batch predictions
  disagreed with row-by-row predictions. The mapping from the
  random-effect model matrix back to a grouping level used
  [`max.col()`](https://rdrr.io/r/base/maxCol.html), which returns a
  column even for an all-zero row (an unobserved interaction level) and
  breaks the resulting tie *at random* — so an unobserved level silently
  borrowed a randomly chosen observed level’s random effect. All-zero
  rows are now detected and routed through the existing new-level path,
  which zeroes that random-effect term’s contribution so the prediction
  falls back to the fixed effects plus any observed higher-level random
  effects, exactly as the out-of-sample warning describes.
  Observed-level predictions are bit-for-bit unchanged.

- **Restored single RNG stream in
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md).**
  The refactor of
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  into component helpers inadvertently had each helper
  ([`simulate_random_effects()`](https://jknowles.github.io/merTools/reference/simulate_random_effects.md),
  [`simulate_fixed_effects()`](https://jknowles.github.io/merTools/reference/simulate_fixed_effects.md))
  call `set.seed(seed)` internally when invoked from the main function,
  resetting the RNG stream mid-call. This produced numerically different
  output for any user-supplied seed compared to CRAN releases. Fixed by
  passing `seed = NULL` from
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  to the helpers; the outer `set.seed(seed)` now pins a single,
  sequential stream (sigma → random effects → fixed effects → residuals)
  exactly as before the refactor. Verified against `origin/master` with
  a 16-case harness: 14 of 16 cases are now bit-for-bit identical to
  master, and the remaining 2 differ only by the intentional GLMM
  binomial-residual simulation fix below.

- Fixed parse error caused by debug statements inserted mid-function
  call in
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  (now correctly handles the function call structure)

- Fixed GLMM residual variance simulation to properly return NULL for
  all GLMM/NLMM models (no Gaussian residual variance exists for
  discrete distributions)

- Added `!is.null(sigma_vec)` checks in
  [`combine_components()`](https://jknowles.github.io/merTools/reference/combine_components.md)
  to handle cases where GLMMs don’t have Gaussian residual variance to
  add (prevents errors with `sigma_vec = NULL`)

- Removed `seed` parameter from
  [`simulate_residual_variance()`](https://jknowles.github.io/merTools/reference/simulate_residual_variance.md)
  (seeds are now handled at the
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  level for reproducibility)

- Updated test expectations to reflect that GLMMs return NULL from
  [`simulate_residual_variance()`](https://jknowles.github.io/merTools/reference/simulate_residual_variance.md)

- Replaced bare `subbars()` with
  [`reformulas::subbars()`](https://rdrr.io/pkg/reformulas/man/subbars.html)
  in the random-effect prediction path to resolve the deprecation
  warning from lme4, which has migrated `subbars` to the `reformulas`
  package.

- **Fixed
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  on models with multiple random-effect term blocks per grouping factor
  ([\#118](https://github.com/jknowles/merTools/issues/118)).** Models
  using the double-bar syntax (`(x + y || g)`), explicit splits
  (`(1|g) + (0 + x|g)`), or mixed correlated + uncorrelated specs
  previously failed with
  `Error in dimnames(reMatrix) <- *vtmp* : 'dimnames' applied to non-array`
  because `lme4::ranef(..., condVar = TRUE)` returns `postVar` as a list
  of per-block arrays in those cases, and the level-filtering code
  assumed a 3-D array.
  [`simulate_random_effects()`](https://jknowles.github.io/merTools/reference/simulate_random_effects.md)
  now normalizes the list to a single block-diagonal array (zero
  off-diagonals between uncorrelated blocks, preserving full covariance
  within correlated blocks) before indexing, so the
  [`mvtnorm::rmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)
  path sees the mathematically correct joint posterior covariance.
  Correlated-only models are unaffected (guarded by
  [`is.list()`](https://rdrr.io/r/base/list.html) check).

- **Fixed
  [`averageObs()`](https://jknowles.github.io/merTools/reference/averageObs.md)
  /
  [`findFormFuns()`](https://jknowles.github.io/merTools/reference/findFormFuns.md)
  on matrix-LHS models
  ([\#83](https://github.com/jknowles/merTools/issues/83)).**
  `averageObs(gm1)` previously errored on two-column binomial GLMMs such
  as `glmer(cbind(successes, failures) ~ ..., family = binomial)`
  because
  [`collapseFrame()`](https://jknowles.github.io/merTools/reference/collapseFrame.md)
  attempted to take the mean of the matrix response column, and a latent
  bug in the weights-selection path tried to index a `(weights)` column
  that does not exist in the model frame for cbind specifications.
  Matrix response columns are now detected and dropped before averaging.
  **Behavior change:** for matrix-LHS models the returned frame no
  longer contains the response column; for scalar-LHS models the
  response is still included as before. Callers that key off column
  count or column names from
  [`averageObs()`](https://jknowles.github.io/merTools/reference/averageObs.md)
  should treat matrix-LHS output as predictors-only. The output remains
  valid `newdata` for
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md),
  which ignore the response column.

### Technical Improvements

- The residual variance logic now correctly distinguishes between:
  - **LMMs**: Gaussian residual variance via `rnorm(N, yhat, sigma)`
    from gamma-distributed sigma
  - **GLMMs**: Conditional distribution simulation via
    [`simulate_glmm_response()`](https://jknowles.github.io/merTools/reference/simulate_glmm_response.md)
    (binomial/poisson/gamma)
- When `include.resid.var = TRUE` and GLMM with `type = "probability"`:
  Simulates from conditional distribution (theoretically correct for
  discrete distributions)
- When `include.resid.var = TRUE` and GLMM with
  `type = "linear.prediction"`: Returns linear predictor without
  Gaussian noise (correct behavior - GLMMs don’t have additive Gaussian
  noise)

### Test Infrastructure

- **Added `tests/comparisons/predictInterval-regression.R`**, a
  standalone cross-version numeric regression harness. It pins a
  canonical set of LMM and GLMM inputs — covering `which`, `level`,
  `stat`, `ignore.fixed.terms`, `fix.intercept.variance`, and
  single-row-newdata cases — and serializes
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  output to an RDS bundle so two package versions can be compared
  bit-for-bit via a `diff` subcommand. Invoke it whenever touching
  simulation internals to confirm the change does not silently alter
  user-facing numeric output. See the README for a worked
  `git worktree`-based workflow.
- Pinned RNG version and algorithm across R releases via a new
  `tests/testthat/helper-seed.R` that sets `RNGversion("4.1.0")` and an
  explicit [`RNGkind()`](https://rdrr.io/r/base/Random.html). This
  prevents silent stream differences across R-oldrel / R-release /
  R-devel on CI.
- Disabled testthat parallel execution
  (`Config/testthat/parallel: false`). File-level
  [`set.seed()`](https://rdrr.io/r/base/Random.html) behaves
  unpredictably under parallel workers with separate RNG state, which
  was the root cause of several of the intermittent CI failures observed
  during the 0.9.0 development cycle.
- Unified all test seeds to `11213` for consistency, preserving a single
  differing seed in two tests that explicitly assert that different
  seeds produce different results.
- Refactored the
  [`thetaExtract()`](https://jknowles.github.io/merTools/reference/thetaExtract.md)
  test in `test-helpers.R` from a brittle numeric-equality check against
  a value calibrated to an older seed, into behavioral assertions (type,
  length, bounds).

### CI

- Bumped `actions/checkout` to `@v5` and set
  `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24: true` on both workflows to
  address the GitHub Actions Node.js 20 deprecation (scheduled for
  2026-09-16).

## merTools 0.6.5

### Code Architecture Improvements

- Refactored
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  into modular component functions for improved maintainability and
  testability. The main function now orchestrates five internal helper
  functions:
  - [`simulate_residual_variance()`](https://jknowles.github.io/merTools/reference/simulate_residual_variance.md) -
    Draws residual standard deviation samples from the posterior
  - [`simulate_fixed_effects()`](https://jknowles.github.io/merTools/reference/simulate_fixed_effects.md) -
    Simulates fixed effect predictions with proper variance-covariance
    handling
  - [`simulate_random_effects()`](https://jknowles.github.io/merTools/reference/simulate_random_effects.md) -
    Simulates random effect contributions for all grouping factors
  - [`combine_components()`](https://jknowles.github.io/merTools/reference/combine_components.md) -
    Combines fixed, random, and residual variance components
  - [`summarise_predictions()`](https://jknowles.github.io/merTools/reference/summarise_predictions.md) -
    Computes prediction intervals from simulation results
- This refactoring reduces the main
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  function from ~520 lines to ~180 lines while preserving complete
  backward compatibility
- Added comprehensive unit tests for all helper functions (43 new tests)
- All existing tests pass without modification, ensuring numeric
  accuracy is preserved

### Benefits

- **Easier maintenance**: Each component function has a single
  responsibility and can be updated independently
- **Better testability**: Individual simulation components can now be
  unit tested in isolation
- **Improved readability**: The main
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  function now clearly shows the high-level algorithm flow
- **Foundation for future enhancements**: The modular architecture makes
  it easier to add new features or optimization strategies

### Notes

- The helper functions are internal and not exported; the public API
  remains unchanged
- Parallelization support is preserved in
  [`simulate_random_effects()`](https://jknowles.github.io/merTools/reference/simulate_random_effects.md)
- Seed reproducibility is maintained by setting the random seed once at
  the start of
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)

## merTools 0.6.4

CRAN release: 2026-01-23

- Maintenance release to merge
  [@DavisVaughan](https://github.com/DavisVaughan) changes to
  accommodate upstream changes in `vctrs` package impacting
  [`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
  usage in `REsim`
  ([\#133](https://github.com/jknowles/merTools/issues/133))

## merTools 0.6.3

CRAN release: 2025-09-05

- Maintenance release to fix crossreference issues with function
  documentation

## merTools 0.6.2

CRAN release: 2024-02-08

- Maintenance release to fix minor issues with function documentation
- Fix [\#130](https://github.com/jknowles/merTools/issues/130) by
  avoiding conflict with `vcov` in the `merDeriv` package
- Upgrade package test infrastructure to 3e testthat specification

## merTools 0.6.1

CRAN release: 2023-03-20

- Maintenance release to keep package listed on CRAN
- Fix a small bug where parallel code path is run twice
  ([\#126](https://github.com/jknowles/merTools/issues/126))
- Update plotting functions to avoid deprecated
  [`aes_string()`](https://ggplot2.tidyverse.org/reference/aes_.html)
  calls ([\#127](https://github.com/jknowles/merTools/issues/127))
- Fix ([\#115](https://github.com/jknowles/merTools/issues/115)) in
  description
- Speed up PI using [@bbolker](https://github.com/bbolker) pull request
  ([\#120](https://github.com/jknowles/merTools/issues/120))
- Updated package maintainer contact information

## merTools 0.5.2

CRAN release: 2020-06-23

- Streamline vignette building to be precompiled and move tests to limit
  burden on CRAN check
- Switch dependency from `broom` to `broom.mixed` because of upstream
  package reorganization

## merTools 0.5.1

### Bug fixes

- Fixed an issue where `averageObs` could not be calculated when model
  weights were specified in the original model (closes
  [\#110](https://github.com/jknowles/merTools/issues/110))

## merTools 0.5.0

CRAN release: 2019-05-13

### New Features

- `subBoot` now works with `glmerMod` objects as well
- `reMargins` a new function that allows the user to marginalize the
  prediction over breaks in the distribution of random effect
  distributions, see `?reMargins` and the new `reMargins` vignette
  (closes [\#73](https://github.com/jknowles/merTools/issues/73))

### Bug fixes

- Fixed an issue where known convergence errors were issuing warnings
  and causing the test suite to not work
- Fixed an issue where models with a random slope, no intercept, and no
  fixed term were unable to be predicted
  ([\#101](https://github.com/jknowles/merTools/issues/101))
- Fixed an issue with shinyMer not working with substantive fixed
  effects ([\#93](https://github.com/jknowles/merTools/issues/93))

## merTools 0.4.2

### New Features

- Parallel fitting of `merModLists` is now supported using the
  `future.apply` package and the `future_lapply` functions, optionally
- Reduced package installation surface by eliminating unnecessary
  packages in the `Suggests` field

### Bug fixes

- Fixed a bug ([\#94](https://github.com/jknowles/merTools/issues/94))
  where
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  would return a data.frame of the wrong dimensions when predicting a
  single row of observations for a `glm`
- Fixed a bug ([\#96](https://github.com/jknowles/merTools/issues/96))
  related to `rstanarm` dependencies in the package vignette
- Switched from `dontrun` to `donttest` for long-running examples (CRAN
  compliance)
- Fixed and made more clear the generics applying to `merModList`
  objects ([\#92](https://github.com/jknowles/merTools/issues/92))

## merTools 0.4.1

CRAN release: 2018-06-05

### New Features

- Standard errors reported by `merModList` functions now apply the Rubin
  correction for multiple imputation

### Bug fixes

- Contribution by Alex Whitworth
  ([@alexWhitworth](https://github.com/alexWhitworth)) adding error
  checking to plotting functions
- The vignettes have been shortened and unit tests reorganized to
  facilitate Travis-CI builds and reduce CRAN build burden

## merTools 0.4.0

### New Features

- Added vignette on using multilevel models with multiply imputed data
- Added `fixef` and `ranef` generics for `merModList` objects
- Added `fastdisp` generic for `merModList`
- Added `summary` generic for `merModList`
- Added `print` generic for `merModList`
- Documented all generics for `merModList` including examples and a new
  imputation vignette
- Added `modelInfo` generic for `merMod` objects that provides simple
  summary stats about a whole model

### Bug Fixes

- Fix bug that returned NaN for `std.error` of a multiply imputed
  `merModList` when calling `modelRandEffStats`
- Fixed bug in `REimpact` where some column names in `newdata` would
  prevent the prediction intervals from being computed correctly. Users
  will now be warned.
- Fixed bug in `wiggle` where documentation incorrectly stated the
  arguments to the function and the documentation did not describe
  function correctly

## merTools 0.3.1

- Update the `readme.rmd` to package graphics with the R package, per
  CRAN

## merTools 0.3.0

CRAN release: 2016-12-12

- Improve handling of formulas. If the original `merMod` has functions
  specified in the formula, the `draw` and `wiggle` functions will check
  for this and attempt to respect these variable transformations. Where
  this is not possible a warning will be issued. Most common
  transformations are respected as long as the the original variable is
  passed untransformed to the model.
- Change the calculations of the residual variance. Previously residual
  variance was used to inflate both the variance around the fixed
  parameters and around the predicted values themselves. This was
  incorrect and resulted in overly conservative estimates. Now the
  residual variance is appropriately only used around the final
  predictions
- Rebuilt the readme.md to include new information about new features
- New option for `predictInterval` that allows the user to return the
  full interval, the fixed component, the random component, or the fixed
  and each random component separately for each observation
- Fixed a bug with slope+intercept random terms that caused a
  miscalculation of the random component
- Add comparison to `rstanarm` to the Vignette
- Make `expectedRank` output more `tidy` like and allow function to
  calculate expected rank for all terms at once
  - Note, this breaks the API by changing the names of the columns in
    the output of this function
- Remove tests that test for timing to avoid issues with R-devel JIT
  compiler
- Remove `plyr` and replace with `dplyr`
- Fix issue [\#62](https://github.com/jknowles/merTools/issues/62)
  `varList` will now throw an error if `==` is used instead of `=`
- Fix issue [\#54](https://github.com/jknowles/merTools/issues/54)
  `predictInterval` did not included random effects in calculations when
  `newdata` had more than 1000 rows and/or user specified
  `parallel=TRUE`. Note: fix was to disable the `.paropts` option for
  `predictInterval` … user can still specify for *temporary* backward
  compatibility but this should be either removed or fixed in the
  permanent solution.
- Fix issue [\#53](https://github.com/jknowles/merTools/issues/53) about
  problems with `predictInterval` when only specific levels of a
  grouping factor are in `newdata` with the colon specification of
  interactions
- Fix issue [\#52](https://github.com/jknowles/merTools/issues/52) ICC
  wrong calculations … we just needed to square the standard deviations
  that we pulled

## merTools 0.2.1

CRAN release: 2016-03-30

- Fix dependency on `lme4` to ensure compatibility with latest changes.

## merTools 0.2

### Bug fixes

- Coerce `dplyr` `tbl` and `tbl_df` objects to data.frames when they are
  passed to `predictInterval` and issue a warning
- Try to coerce other data types passed to `newdata` in
  `predictInterval` before failing if coercion is unsuccessful
- Numeric stabilization of unit tests by including seed values for
  random tests
- Fix handling of models with nested random effect terms (GitHub
  [\#47](https://github.com/jknowles/merTools/issues/47))
- Fix vignette images

### New Functionality

- Substantial performance enhancement for `predictInterval` which
  includes better handling of large numbers of parameters and
  simulations, performance tweaks for added speed (~10x), and parallel
  backend support (currently not optimized)
- Add support for `probit` models and limited support for other `glmm`
  link functions, with warning (still do not know how to handle sigma
  parameter for these)
- Add ability for user-specified seed for reproducibility
- Add support for `blmer` objects from the `blme` package
- Add a `merModList` object for lists of `merMod` objects fitted to
  subsets of a dataset, useful for imputation or for working with
  extremely large datasets
- Add a `print` method for `merModList` to mimic output of
  `summary.merMod`
- Add a `VarCorr` method for `merModList`
- Add new package data to demonstrate replication from selected
  published texts on multilevel modeling using different software (1982
  High School and Beyond Survey data)

### Other changes

- Changed the default `n.sims` for the `predictInterval` function from
  100 to 1,000 to give better coverage and reflect performance increase
- Changed the default for `level` in `predictInterval` to be 0.8 instead
  of 0.95 to reflect that 0.95 prediction intervals are more
  conservative than most users need

### Future changes

- For the next release (1.0) we are considering a permanent switch to
  C++ RMVN sampler courtesy of Giri Gopalan ’s excellent FastGP

## merTools 0.1

- Initial release

### New Functions

- Provides `predictInterval` to allow prediction intervals from `glmer`
  and `lmer` objects
- Provides `FEsim` and `REsim` to extract distributions of model
  parameters
- Shows `shinyMer` an interactive `shiny` application for exploring
  `lmer` and `glmer` models
- Provides `expectedRank` function to interpret the ordering of effects
- Provides `REimpact` to simulate the impact of grouping factors on the
  outcome
- Provides `draw` function to allow user to explore a specific
  observation
- Provides `wiggle` function for user to build a simulated set of
  counterfactual cases to explore
