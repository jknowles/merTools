# merTools 1.0 Release TODO

*Created: 2026-05-23. Tracks remaining work from
`plans/1.0-release-readiness.md`.* *Current version: 0.9.0. Target:
1.0.0 on CRAN.*

------------------------------------------------------------------------

## Phase 1 — Housekeeping sprint (2026-05-23)

Close GitHub issues already fixed in code: \#118, \#121 (closed); \#83,
\#101, \#134 (already closed)

Run `urlchecker::url_check()` — fixed JSTOR URL in `expectedRank` (→
`\doi{10.2307/1164724}`)

Run `spelling::spell_check_package()` — fixed 5 typos (accommodate,
identifiable, substantive, programmatically, co-occurrences); added
`inst/WORDLIST` for technical terms

Run `covr::package_coverage()` — **baseline: 46.85%** (below 80%
target). 8 files at 0% — likely test-skipping issue, needs investigation

Decide `residDF.merMod` scale-parameter — **keep current behavior**
(scale param counted). Matches Gelman & Hill (2007). TODO removed,
rationale documented in code.

Decide roxygen markdown mode — **stay LaTeX for 1.0**. ~83 `\code{}`
instances; conversion deferred to 1.1

Fix test hygiene: removed `rm(list = ls())` (3 occurrences), replaced
`skip_on_travis()` → `skip_on_ci()` (2), replaced `sink("NUL")` →
`sink(nullfile())` (1), deleted NUL artifact

Triage \#136, \#137 — docs/citation improvements, folded into Phase 4

Regenerate Rd files (expectedRank.Rd, REmargins.Rd,
combine_postvar_blocks.Rd)

### Coverage investigation needed

The 46.85% baseline has 8 source files at 0% despite tests existing for
several (merFastDisplay, merList, merPlots, merSubstEff, REmargins).
Possible causes: - Split test runners (`filter = "^[a-m]"` / `"^[m-z]"`)
may miss files - Tests may fail silently in covr’s environment - This
should be investigated before Phase 3 test additions

## Phase 2 — Issue \#124 investigation

Reproduce \#124: non-reproducible `predictInterval` with mixed
observed/unobserved levels

Determine root cause (mkNewReTrms, reMatrix path, or FE simulation)

Fix or document as expected behavior in
[`?predictInterval`](reference/predictInterval.md)

Add regression test for the chosen behavior

## Phase 3 — API cleanup

Import hygiene: convert `import(dplyr)` to specific `importFrom` calls

Import hygiene: convert `import(ggplot2)` to specific `importFrom` calls

Import hygiene: convert `import(lme4)` to specific `importFrom` calls

Import hygiene: convert `import(arm)` to specific `importFrom` calls
(already partially done)

Add `lifecycle` to Imports

Annotate all 40 exports with lifecycle badges
(stable/experimental/deprecated)

Add basic tests for untested exports: modelInfo, modelFixedEff,
modelRandEffStats

Add basic tests for untested exports: REsdExtract, REcorrExtract,
randomObs

Add basic tests for untested exports: bglmerModList, fixef.merModList,
ranef.merModList

Document [`REmargins()`](reference/REmargins.md) tie-breaking behavior;
add performance bounds note

Investigate and fix test coverage gaps (8 files at 0%)

## Phase 4 — Documentation pass

Refresh vignettes: update dates, verify code runs, re-knit via .Rmd.orig

Fix bootMer graphic alignment in Using_predictInterval vignette (#116)

Update README.Rmd: add 0.9.0 and 1.0 entries to “Recent Updates”

Re-knit README.md from README.Rmd

Write NEWS.md 1.0 section: version jump rationale, API stability pledge,
breaking changes

Decide on shinyMer landing tab (#78) — add Model Summary tab or defer to
post-1.0

Address \#136 (fresh examples from easystats/modelbased) and \#137
(citation improvements)

## Phase 5 — Final release checks

Run `revdepcheck::revdep_check()`

Run `devtools::check(remote = TRUE, manual = TRUE)`

Run `R CMD check --as-cran` on the built tarball

Bump version to 1.0.0 in DESCRIPTION

Update cran-comments.md for 1.0 submission

Tag v1.0.0 and create GitHub release

Submit to CRAN

## Deferred to post-1.0

- \#84, \#85 — REimpact plot improvements
- \#32 — shinyMer user-controlled case selection
- Cross-platform Layer 2 snapshots (macOS-ARM64, Windows)
- foreach → future.apply migration
- cli/rlang error modernization
- pkgdown site verification
- roxygen markdown mode conversion
- REmargins performance benchmarking and future backend
