## Submission summary

This is merTools 1.0.0, the 1.0 long-term-support release. It resolves the
remaining open issues, fixes a correctness bug in `predictInterval()` for
nested/interaction random effects, repairs and extends the `shinyMer()`
explorer, adds a new `plotREimpact()` visualization, refreshes the
documentation, and tidies the test suite for low-maintenance, long-term use.
The public API is unchanged and existing results are preserved (the new
`new.levels` behavior is opt-in; the default reproduces prior output).

## Test environments

* local: Ubuntu 26.04 LTS, R 4.6.0
* GitHub Actions:
    - macOS-latest,   R release
    - windows-latest, R release
    - ubuntu-latest,  R devel
    - ubuntu-latest,  R release
    - ubuntu-latest,  R oldrel-1
* win-builder: R devel and R release

## R CMD check results

0 errors | 0 warnings | 0 notes

## revdepcheck results

We checked 3 reverse dependencies, comparing R CMD check results across CRAN and
dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
