## Test environments
* local Windows 11 install, R 4.3.1
* r-hub (R-release macOS 10.11, R-devel)
* win-builder (R-release, R-devel)

## R CMD check results

* Updated \usage sections
* Used `utils::globalVariables(c(".shinyMerPar", "sig", "sigma"))` to fix notes
  about unexported objects.

* Used `utils::globalVariables(c("Lind", "group", "est", "mean_est", "est_ss", 
                                "within_var", "between_var", "statistic"))` 
    to fix notes about unexported objects in `modelFixedEff` and 
    `fastdisp.merList` functions
    

## Downstream dependencies
There are currently two downstream dependencies. No issues with either were 
found.
