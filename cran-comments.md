## Test environments

* local Ubuntu linux install, R 4.5.2
* GitHub Actions
    - {os: macos-latest,   r: 'release'}
    - {os: windows-latest, r: 'release'}
    - {os: ubuntu-latest,   r: 'devel', http-user-agent: 'release'}
    - {os: ubuntu-latest,   r: 'release'}
    - {os: ubuntu-latest,   r: 'oldrel-1'}

## R CMD check results

* Fixed all issues with upstream dependencies caused by `dplyr::bind_rows()`
* Used `utils::globalVariables(c(".shinyMerPar", "sig", "sigma"))` to fix notes
  about unexported objects.
* Used `utils::globalVariables(c("Lind", "group", "est", "mean_est", "est_ss",
                                "within_var", "between_var", "statistic"))`
    to fix notes about unexported objects in `modelFixedEff` and
    `fastdisp.merList` functions

## Reverse depencency check results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
