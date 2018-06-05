## Test environments
* local Windows 10 install, R 3.5.0
* ubuntu 14.05 (on travis-ci), R devel, R-Release
* win-builder (devel and release)

## R CMD check results
* Used utils::globalVariables(c(".shinyMerPar", "sig", "sigma")) to fix notes
about unexported objects.
* Used utils::globalVariables(c("Lind", "group", "est", "mean_est", "est_ss", 
                                "within_var", "between_var", "statistic")) 
to fix notes about unexported objects in `modelFixedEff` and `fastdisp.merList` 
functions
* The Vignettes have been shortened in an attempt to reduce build times

## Downstream dependencies
There are currently two downstream dependencies. No issues with either were found.
