# NEWS

## merTools 0.2

### Bug fixes

- Coerce `dplyr` `tbl` and `tbl_df` objects to data.frames when they are passed 
to `predictInterval` and issue a warning
- Try to coerce other data types passed to `newdata` in `predictInterval` before 
failing if coercion is unsuccessful

### New Functionality

- Performance enhancement for `predictInterval` which includes better 
handling of large numbers of parameters and simulations, performance 
tweaks for added speed (~25%), and parallel backend support
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
