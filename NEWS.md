# NEWS

## merTools 1.1

### Bug fixes

- Coerce `dplyr` `tbl` and `tbl_df` objects to data.frames when they are passed 
to `predictInterval` and issue a warning
- Try to coerce other data types passed to `newdata` in `predictInterval` before 
failing if coercion is unsuccessful

### New Functionality

- Add support for `probit` models and limited support for other `glmm` link functions, with warning
- Add ability for user-specified seed for reproducibility
- Add support for `blmer` objects from the `blme` package

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
