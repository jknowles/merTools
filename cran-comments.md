## Test environments
* local Windows 10 install, R 3.5.3
* r-hub (R-release macOS 10.11)
* win-builder (R-release, R-devel)

## R CMD check results

* Used `utils::globalVariables(c(".shinyMerPar", "sig", "sigma"))` to fix notes
  about unexported objects.

* Used `utils::globalVariables(c("Lind", "group", "est", "mean_est", "est_ss", 
                                "within_var", "between_var", "statistic"))` 
    to fix notes about unexported objects in `modelFixedEff` and 
    `fastdisp.merList` functions
    
* Note about possible invalid URL is a false positive: http://www.jstor.org/stable/1164724
    is a valid URL

## Downstream dependencies
There are currently two downstream dependencies. No issues with either were 
found.
