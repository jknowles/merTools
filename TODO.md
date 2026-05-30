# merTools 1.0 Release TODO

*Updated: 2026-05-29. Version: **1.0.0** (in DESCRIPTION). Target: 1.0.0
on CRAN.*

The 1.0 readiness work (Phases 1–4 below) is **complete**. What remains
is the CRAN submission itself and the release announcement — see the
next section.

------------------------------------------------------------------------

## 🚀 Tomorrow — CRAN Submission Day (2026-05-30)

### A. Pre-flight (make the source release final)

Confirm `DESCRIPTION` `Version: 1.0.0` and `Date` (add/refresh `Date:`
if used).

**Rebuild `README.md` with `devtools::build_readme()`** (NOT bare
`knit()` — that reintroduces the YAML front matter we just stripped) so
badges, results, and figures are clean and current.

Final `devtools::document()`; confirm `RoxygenNote: 7.3.3` unchanged.

Re-read `NEWS.md` 1.0.0 section top-to-bottom; confirm every bullet is
accurate and the LTS / maintenance-mode framing reads well.

Verify `inst/CITATION` renders with version 1.0.0
([`utils::readCitationFile`](https://rdrr.io/r/utils/citation.html)).

`git status` clean; everything committed and pushed to `main`.

### B. Checks (all should be clean before submitting)

`devtools::check()` locally → **0 errors / 0 warnings / 0 notes**.

`R CMD check --as-cran` on the built tarball **with vignettes built**
(confirms the precompiled brms vignette renders without brms installed).

`devtools::spell_check()` (WORDLIST exists — add any new technical
terms, e.g. brms, posterior, GLMM, predictInterval).

`urlchecker::url_check()` (we added a `URL:` field — confirm it
resolves; the pkgdown site is live at
<https://jknowles.github.io/merTools/>).

win-builder: `devtools::check_win_devel()`, `check_win_release()`,
`check_win_oldrelease()` — wait for the 3 email results.

R-hub: `rhub::rhub_check()` (or `rhub::check_for_cran()`) across the
standard CRAN platforms.

Reverse dependencies: `install.packages("revdepcheck")` then
`revdepcheck::revdep_check(num_workers = 4)`. There are **3 revdeps**
(last submission saw 0 new problems) — expect 0 new problems.

### C. cran-comments.md

Update `cran-comments.md` for 1.0.0: - Test environments (local + GitHub
Actions matrix + win-builder/R-hub). - One line: “This is the 1.0
long-term-support release.” - R CMD check results: 0 errors / 0 warnings
/ 0 notes. - Reverse dependency results: “Checked 3 reverse
dependencies; 0 new problems.”

### D. Submit

`devtools::submit_cran()` (writes `CRAN-SUBMISSION`) **or** upload the
tarball at <https://cran.r-project.org/submit.html>.

Confirm the submission via the auto-generated email.

Watch for CRAN’s automated check results / any maintainer feedback;
address promptly and resubmit if needed.

### E. After acceptance

Tag `v1.0.0` and create a GitHub release with the NEWS highlights.

Confirm the CRAN version badge and pkgdown site reflect 1.0.0.

Publish the blog post (section below) and cross-post.

------------------------------------------------------------------------

## 📣 Blog post — “merTools 1.0” announcement

Draft, review, and publish on the Civilytics / jaredknowles.com blog.

Cross-post: R-bloggers (via the blog’s R feed/tag), Mastodon/Bluesky,
rOpenSci/Rweekly if appropriate.

**Starter outline (talking points):**

1.  **Hook** — merTools has been on CRAN since 2015 (~457K downloads,
    ~8K/month). 1.0 marks it *feature complete* and the start of a
    long-term-support phase.
2.  **What it does** — fast, simulation-based prediction intervals for
    `lme4` mixed models
    ([`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)),
    when `bootMer()` / full MCMC is impractical; plus tools to extract,
    simulate, and visualize effects.
3.  **What’s new in 1.0** (keep tight, link to NEWS):
    - Correctness fix: reproducible
      [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
      for nested/interaction random effects (#124).
    - New `new.levels = "draw"` — sample unobserved groups’ effects from
      the RE covariance (matches brms `allow_new_levels`).
    - New
      [`plotREimpact()`](https://jknowles.github.io/merTools/reference/plotREimpact.md);
      [`plotFEsim()`](https://jknowles.github.io/merTools/reference/plotFEsim.md)
      highlights significant terms.
    - [`shinyMer()`](https://jknowles.github.io/merTools/reference/shinyMer.md)
      revived and extended (model-summary tab, subset draws).
    - Refreshed docs: a new *Contextual Effects* vignette and proper
      citation.
4.  **The headline** — *validated against brms.* predictInterval
    reproduces the full-Bayesian posterior prediction intervals almost
    exactly (point estimates `cor 0.9998`, matching 90% intervals and
    coverage) at ~1/400th the cost. Link/screenshot the new “Validating
    predictInterval() against brms” vignette and the calibration figure.
5.  **Why LTS / maintenance mode** — stable API, a decade of use; future
    work is bug fixes, CRAN/dependency compatibility, and docs, not
    churn.
6.  **Call to action** — `install.packages("merTools")`, read the site
    (<https://jknowles.github.io/merTools/>), cite it
    (`citation("merTools")`), file issues. Thank contributors (Carl
    Frederick, Alex Whitworth, bbolker, DavisVaughan, and reporters like
    dotPiano for \#124).

------------------------------------------------------------------------

## Optional before 1.0 (defer to 1.0.x if time-constrained — none are CRAN blockers)

These are Phase 3 polish items; the package passes `--as-cran` without
them.

Import hygiene: `import(pkg)` → specific `importFrom` (dplyr, ggplot2,
lme4, arm).

Add `lifecycle` and badge the ~40 exports (stable/experimental).

Add basic tests for currently-untested exports (modelInfo,
modelFixedEff, modelRandEffStats, REsdExtract, REcorrExtract, randomObs,
merModList generics).

Investigate the coverage gap (46.85% baseline; 8 files reported 0% —
likely a covr/test-runner artifact, not missing tests).

------------------------------------------------------------------------

## Completed in the 1.0 readiness work (for the record)

- **Phase 1 — Housekeeping** (2026-05-23): urlchecker, spell check +
  WORDLIST, test hygiene, residDF/roxygen decisions, Rd regeneration.
- **Phase 2 — \#124**: fixed non-reproducible
  [`predictInterval()`](https://jknowles.github.io/merTools/reference/predictInterval.md)
  for nested / interaction random effects; regression tests added.
- **Phase 4 — Documentation**: vignettes refreshed and re-knit; bootMer
  graphic fixed (#116); README updated (Recent Updates + plotREimpact
  section); NEWS 1.0 section with LTS framing; citation + influences
  (#137); contextual-effects vignette (#136); shinyMer model-summary tab
  (#78).
- **Pulled forward from “deferred”**:
  [`plotREimpact()`](https://jknowles.github.io/merTools/reference/plotREimpact.md)
  and
  [`plotFEsim()`](https://jknowles.github.io/merTools/reference/plotFEsim.md)
  polish (#84, \#85); shinyMer subset/case selection (#32); pkgdown site
  verified and live; brms validation vignette + `new.levels = "draw"`
  feature.

## Still deferred to post-1.0

- foreach → future.apply migration; cli/rlang error modernization.
- roxygen markdown-mode conversion (~83 `\code{}` instances).
- Cross-platform Layer-2 snapshots (macOS-ARM64, Windows).
- `REmargins` performance benchmarking / future backend.
- A bayesplot/ggridges density display for FE/RE plots (noted in \#85).
