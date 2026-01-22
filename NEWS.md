# NEWS

## merTools 0.6.4

- Maintenance release to merge @DavisVaughan changes to accomodate upstream changes in `vctrs` package
- Fixes changes in `findbars` and `mkReTrms` upstream


## merTools 0.6.3

- Maintenance release to fix crossreference issues with function documentation

## merTools 0.6.2

- Maintenance release to fix minor issues with function documentation
- Fix #130 by avoiding conflict with `vcov` in the `merDeriv` package
- Upgrade package test infrastructure to 3e testthat specification

## merTools 0.6.1

- Maintenance release to keep package listed on CRAN
- Fix a small bug where parallel code path is run twice (#126)
- Update plotting functions to avoid deprecated `aes_string()` calls (#127)
- Fix (#115) in description
- Speed up PI using @bbolker pull request (#120)
- Updated package maintainer contact information

## merTools 0.5.2

- Streamline vignette building to be precompiled and move tests to limit burden on CRAN check
- Switch dependency from `broom` to `broom.mixed` because of upstream package reorganization

## merTools 0.5.1

### Bug fixes

- Fixed an issue where `averageObs` could not be calculated when model weights were specified in the
original model (closes #110)

## merTools 0.5.0

### New Features

- `subBoot` now works with `glmerMod` objects as well
- `reMargins` a new function that allows the user to marginalize the prediction over breaks in the
distribution of random effect distributions, see `?reMargins` and the new `reMargins` vignette (closes #73)

### Bug fixes

- Fixed an issue where known convergence errors were issuing warnings and causing the test suite
to not work
- Fixed an issue where models with a random slope, no intercept, and no fixed term were unable
to be predicted (#101)
- Fixed an issue with shinyMer not working with substantive fixed effects (#93)


## merTools 0.4.2

### New Features

- Parallel fitting of `merModLists` is now supported using the `future.apply`
package and the `future_lapply` functions, optionally
- Reduced package installation surface by eliminating unnecessary packages
in the `Suggests` field

### Bug fixes

- Fixed a bug (#94) where `predictInterval()` would return a data.frame of the
wrong dimensions when predicting a single row of observations for a `glm`
- Fixed a bug (#96) related to `rstanarm` dependencies in the package vignette
- Switched from `dontrun` to `donttest` for long-running examples (CRAN compliance)
- Fixed and made more clear the generics applying to `merModList` objects (#92)

## merTools 0.4.1

### New Features

- Standard errors reported by `merModList` functions now apply the Rubin
correction for multiple imputation

### Bug fixes
- Contribution by Alex Whitworth (@alexWhitworth) adding error checking to plotting
functions
- The vignettes have been shortened and unit tests reorganized to facilitate
Travis-CI builds and reduce CRAN build burden

## merTools 0.4.0

### New Features
- Added vignette on using multilevel models with multiply imputed data
- Added `fixef` and `ranef` generics for `merModList` objects
- Added `fastdisp` generic for `merModList`
- Added `summary` generic for `merModList`
- Added `print` generic for `merModList`
- Documented all generics for `merModList` including examples and a new
imputation vignette
- Added `modelInfo` generic for `merMod` objects that provides simple summary
stats about a whole model

### Bug Fixes
- Fix bug that returned NaN for `std.error` of a multiply imputed `merModList`
when calling `modelRandEffStats`
- Fixed bug in `REimpact` where some column names in `newdata` would prevent the
prediction intervals from being computed correctly. Users will now be warned.
- Fixed bug in `wiggle` where documentation incorrectly stated the arguments to
the function and the documentation did not describe function correctly

## merTools 0.3.1

- Update the `readme.rmd` to package graphics with the R package, per CRAN

## merTools 0.3.0

- Improve handling of formulas. If the original `merMod` has functions specified
in the formula, the `draw` and `wiggle` functions will check for this and attempt
to respect these variable transformations. Where this is not possible a warning
will be issued. Most common transformations are respected as long as the the
original variable is passed untransformed to the model.
- Change the calculations of the residual variance. Previously residual variance
was used to inflate both the variance around the fixed parameters and around the
predicted values themselves. This was incorrect and resulted in overly conservative
estimates. Now the residual variance is appropriately only used around the
final predictions
- Rebuilt the readme.md to include new information about new features
- New option for `predictInterval` that allows the user to return the full
interval, the fixed component, the random component, or the fixed and each random
component separately for each observation
- Fixed a bug with slope+intercept random terms that caused a miscalculation of
the random component
- Add comparison to `rstanarm` to the Vignette
- Make `expectedRank` output more `tidy` like and allow function to calculate
expected rank for all terms at once
  - Note, this breaks the API by changing the names of the columns in the output
  of this function
- Remove tests that test for timing to avoid issues with R-devel JIT compiler
- Remove `plyr` and replace with `dplyr`
- Fix issue #62 `varList` will now throw an error if `==` is used instead of `=`
- Fix issue #54 `predictInterval` did not included random effects in calculations
  when `newdata` had more than 1000 rows and/or user specified `parallel=TRUE`.
  Note: fix was to disable the `.paropts` option for `predictInterval` ... user
  can still specify for *temporary* backward compatibility but this should be
  either removed or fixed in the permanent solution.
- Fix issue #53 about problems with `predictInterval` when only specific levels
  of a grouping factor are in `newdata` with the colon specification of
  interactions
- Fix issue #52 ICC wrong calculations ... we just needed to square the standard
  deviations that we pulled

## merTools 0.2.1

- Fix dependency on `lme4` to ensure compatibility with latest changes.

## merTools 0.2

### Bug fixes

- Coerce `dplyr` `tbl` and `tbl_df` objects to data.frames when they are passed
to `predictInterval` and issue a warning
- Try to coerce other data types passed to `newdata` in `predictInterval` before
failing if coercion is unsuccessful
- Numeric stabilization of unit tests by including seed values for random tests
- Fix handling of models with nested random effect terms (GitHub #47)
- Fix vignette images

### New Functionality

- Substantial performance enhancement for `predictInterval` which includes better
handling of large numbers of parameters and simulations, performance
tweaks for added speed (~10x), and parallel backend support (currently not optimized)
- Add support for `probit` models and limited support for other `glmm` link functions, with warning (still do not know how to handle sigma parameter
for these)
- Add ability for user-specified seed for reproducibility
- Add support for `blmer` objects from the `blme` package
- Add a `merModList` object for lists of `merMod` objects fitted to subsets
of a dataset, useful for imputation or for working with extremely large datasets
- Add a `print` method for `merModList` to mimic output of `summary.merMod`
- Add a `VarCorr` method for `merModList`
- Add new package data to demonstrate replication from selected published texts
on multilevel modeling using different software (1982 High School and Beyond Survey data)

### Other changes

- Changed the default `n.sims` for the `predictInterval` function from 100 to 1,000
to give better coverage and reflect performance increase
- Changed the default for `level` in `predictInterval` to be 0.8 instead of 0.95
to reflect that 0.95 prediction intervals are more conservative than most users
need

### Future changes
- For the next release (1.0) we are considering a permanent switch to
C++ RMVN sampler courtesy of Giri Gopalan 's excellent FastGP



## merTools 0.1
- Initial release

### New Functions
- Provides `predictInterval` to allow prediction intervals from `glmer` and `lmer`
objects
- Provides `FEsim` and `REsim` to extract distributions of model parameters
- Provides `shinyMer` an interactive `shiny` application for exploring `lmer`
and `glmer` models
- Provides `expectedRank` function to interpret the ordering of effects
- Provides `REimpact` to simulate the impact of grouping factors on the outcome
- Provides `draw` function to allow user to explore a specific observation
- Provides `wiggle` function for user to build a simulated set of counterfactual
cases to explore
