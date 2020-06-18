# Precompiled vignettes that depend on API key
# Must manually move image files from eia/ to eia/vignettes/ after knit

library(knitr)
knit("vignettes/Using_predictInterval.Rmd.orig",
     "vignettes/Using_predictInterval.Rmd")

knit("vignettes/merToolsIntro.Rmd.orig",
     "vignettes/merToolsIntro.Rmd")

knit("vignettes/marginal_effects.Rmd.orig",
     "vignettes/marginal_effects.Rmd")
knit("vignettes/imputation.Rmd.orig",
     "vignettes/imputation.Rmd")

