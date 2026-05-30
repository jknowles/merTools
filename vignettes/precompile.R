# Precompile the vignettes.
#
# These vignettes are *precompiled*: the expensive code (model fitting,
# bootstrapping, Stan sampling) is run here, once, and the resulting .Rmd files
# (with output and figures embedded) are checked in. This keeps `R CMD check`
# and the CRAN build fast and free of heavy Suggested dependencies. Run this
# script from the package root after changing any *.Rmd.orig source, then
# rebuild the package documentation.
#
# The figure paths in each source use a relative fig.path (e.g. "usage-",
# "context-"), so we knit with the working directory set to vignettes/ to keep
# the generated images alongside the .Rmd files.

library(knitr)

local({
  old_wd <- setwd("vignettes")
  on.exit(setwd(old_wd), add = TRUE)

  knit("Using_predictInterval.Rmd.orig", "Using_predictInterval.Rmd")
  knit("merToolsIntro.Rmd.orig",         "merToolsIntro.Rmd")
  knit("marginal_effects.Rmd.orig",      "marginal_effects.Rmd")
  knit("imputation.Rmd.orig",            "imputation.Rmd")
  knit("contextual_effects.Rmd.orig",    "contextual_effects.Rmd")
  # Requires the 'brms' package + a Stan toolchain; only knit when available.
  if (requireNamespace("brms", quietly = TRUE)) {
    knit("brms_validation.Rmd.orig",     "brms_validation.Rmd")
  }
})
