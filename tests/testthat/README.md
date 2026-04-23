# merTools test suite layers

The tests in this directory are organized into three conceptual layers. Each
layer answers a different question. Knowing which layer a test belongs to
determines how brittle it is allowed to be, how fast it must run, and what
platforms it can run on.

When you write a new test, pick the layer that matches the question you are
actually trying to answer. Mixing layers in one test (e.g., bolting a tight
Monte Carlo bias check onto a structural compatibility test) is the most
common cause of flaky tests in this package — don't do it.

---

## Layer 1 — Structural invariants

> *"Is the output still a valid `predictInterval()` result?"*

File: **`test-predictInterval-invariants.R`**

Tests properties that MUST hold for every valid output regardless of seed,
simulation noise, platform, model type, or option combination. Derived from
the function's contract, not from specific numbers.

Examples: output is a `data.frame`, has the expected column names, `upr >=
fit >= lwr` row-wise, `nrow(output) == nrow(newdata)`, no `NA`/`NaN`/`Inf`,
same seed produces identical output, higher `level` gives wider intervals,
`type = "probability"` values are in `[0, 1]`.

**Runs on every CI platform.** No `skip_on_os`. No `skip_on_ci`. Seed- and
platform-independent.

When adding a new invariant: phrase it as "for any valid input, X must hold"
and verify it is derivable from the function's documentation or mathematics.
Do not snapshot numeric output here — that's Layer 2.

---

## Layer 2 — Pinned numeric regression

> *"Did my edit silently change the numbers users see?"*

File: **`test-predictInterval-snapshot.R`**
Committed contract: **`_snaps/predictInterval-snapshot.md`**

Uses `testthat::expect_snapshot_value(style = "json2", tolerance = 1e-6)` to
pin the exact `predictInterval()` output for a canonical set of LMM and
GLMM inputs. Any code change that alters these outputs produces a
reviewable diff in the committed `.md` file, not a random-tolerance failure.

**Runs on Linux only.** macOS-ARM64 and Windows produce bit-level drift from
different BLAS/LAPACK that would pollute the snapshot. Linux is the
canonical platform for the committed contract. Cross-platform equivalence
is covered manually by Layer 3.

When adding a new model to Layer 2: first check that the design matrix is
NOT rank-deficient or near-singular. `mvtnorm::rmvnorm(method = "chol")` on
such a matrix produces output that differs across BLAS implementations by
enough to blow past any reasonable snapshot tolerance, even on the same OS
across R versions. Rank-deficient models belong in Layer 1 and Layer 3 only.

To intentionally accept a changed snapshot:

    testthat::snapshot_accept("predictInterval-snapshot")

Review the resulting `.md` diff carefully before committing.

---

## Layer 3 — Statistical validation and offline regression

> *"Is the algorithm statistically correct (unbiased, well-covered)?
>   And does my refactor preserve numeric behavior across versions?"*

There are two mechanisms for Layer 3 content:

### In-suite, tagged with `skip_on_ci()`

Tests that require large `n.sims`, simulate new data, check coverage or
unbiasedness with tight tolerances, or exercise large datasets. These
still run automatically during local `devtools::test()` but skip on CI,
where they would be slow and flaky by nature.

Examples in this directory:

- `test-predict1.R` — coverage and bias checks on simulated LMM/GLMM data
- `test-predict2.R` — tight-tolerance median-PI accuracy checks
- `test-predict3.R` — rank-deficient model accuracy
- `test-predict4.R` — parallel execution performance/consistency
- `test-predictInterval-helpers.R` (one block) — large `InstEval` fit

When adding a new statistical check, tag it with `skip_on_ci()` and add a
comment explaining what statistical property it is validating.

### Cross-version regression harness (manual)

File: **`tests/comparisons/predictInterval-regression.R`**

A two-mode script (`harness` + `diff`) that pins a canonical set of
`predictInterval()` inputs and serializes their outputs as RDS, enabling
bit-for-bit comparison between two package versions. Run manually when
touching simulation internals to confirm a refactor did not change
numeric behavior for users with a fixed seed.

See the script header or the top-level `README.md` for the full workflow.

---

## Other test files in this directory

The remaining files test **specific model types and behaviors** that are
not numerical regressions or universal invariants:

- `test-predict.R` — edge-case model specifications (nested effects, cross-
  level interactions, no fixed intercept). Structural compatibility only.
- `test-predictInterval-accuracy.R` — reproducibility and dimension checks
  for a range of options.
- `test-predictInterval-helpers.R` — white-box unit tests for internal
  helper functions (`simulate_residual_variance`, `simulate_fixed_effects`,
  etc.).
- `test-seeds.R` — asserts that `seed =` arguments behave correctly across
  all user-facing stochastic functions.
- `test-merData.R`, `test-merExtract.R`, `test-merModList.R`,
  `test-expectedRank.R`, `test-REmargins.R`, `test-substEff.R`, `test-helpers.R`,
  `test-subboot.R` — domain-specific tests for the non-`predictInterval`
  surface of `merTools`.

---

## Conventions

- **Seed**: Use `seed = 11213` (or `set.seed(11213)`) unless a different
  seed is required to demonstrate a specific property (e.g., `test-seeds.R`
  pairs a second seed to assert divergence).
- **`skip_on_cran()`**: Use on any test that takes more than ~1 second.
  CRAN auto-checks are time-sensitive.
- **`skip_on_ci()`**: Use on Layer 3 content per above.
- **`skip_on_os()`**: Use only in `test-predictInterval-snapshot.R` to
  scope Layer 2 to Linux.

## Pinned RNG

`helper-seed.R` sets `RNGversion("4.1.0")` and an explicit `RNGkind()` so
that the same seed produces the same stream across R-oldrel, R-release,
and R-devel on CI. Without this, default RNG tweaks across R versions can
silently change test outputs.
