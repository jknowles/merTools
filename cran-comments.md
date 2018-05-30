## Test environments
* local Windows 10 install, R 3.5.0
* ubuntu 14.05 (on travis-ci), R devel, R-Release
* win-builder (devel and release)

## R CMD check results
There was one NOTE regarding a URL in the `REimpact` function documentation. 
This URL has been checked on all built version of the package and is a valid 
link from the documentation page to the expected target. 

* Used utils::globalVariables(c(".shinyMerPar", "sig", "sigma")) to fix notes
about unexported objects.
* Used utils::globalVariables(c("Lind", "group", "est", "mean_est", "est_ss", 
                                "within_var", "between_var", "statistic")) 
to fix notes about unexported objects in `modelFixedEff` and `fastdisp.merList` 
functions

## Downstream dependencies
There are currently no downstream dependencies. 
