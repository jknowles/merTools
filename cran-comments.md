## Test environments
* local Windows 7 install, R 3.3.2
* ubuntu 12.04 (on travis-ci), R devel, R-Release
* win-builder (devel and release)

## R CMD check results
There was one NOTE regarding a URL in the `REimpact` function documentation. 
This URL has been checked on all built version of the package and is a valid 
link from the documentation page to the expected target. 

* Used utils::globalVariables(c(".shinyMerPar", "sig", "sigma")) to fix notes
about unexported objects.
* Used utils::globalVariables(c("term", "estimate", "std.error")) to fix notes 
about unexported objects in `modelFixedEff` function

## Downstream dependencies
There are currently no downstream dependencies. 
