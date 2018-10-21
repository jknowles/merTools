## Test environments
* local Windows 10 install, R 3.5.1
* win-builder (devel and release)

## R CMD check results
* Warning: parse error in file 'merTools-Ex.R':
  1: unexpected symbol
  119: cleanEx()
  120: nameEx 

This appears to be a false positive warning - all example code runs outside of 
the R CMD CHECK process and the unexpected symbol is related to the parsing as 
part of the R project.

* Used `utils::globalVariables(c(".shinyMerPar", "sig", "sigma"))` to fix notes
  about unexported objects.

* Used `utils::globalVariables(c("Lind", "group", "est", "mean_est", "est_ss", 
                                "within_var", "between_var", "statistic"))` 
    to fix notes about unexported objects in `modelFixedEff` and 
    `fastdisp.merList` functions

## Downstream dependencies
There are currently two downstream dependencies. No issues with either were 
found.
