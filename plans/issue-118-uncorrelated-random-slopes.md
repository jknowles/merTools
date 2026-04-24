# Issue #118 — `predictInterval()` fails on uncorrelated random slopes

## Summary

`predictInterval()` throws
`Error in dimnames(reMatrix) <- *vtmp* : 'dimnames' applied to non-array`
whenever a model has multiple random-effect terms for the same grouping
factor. This occurs with the double-bar syntax (`(x + y || g)`), with
explicit splits (`(1|g) + (0+x|g)`), and with mixed correlated+uncorrelated
specifications (`(Days|Subject) + (0+z|Subject)`).

## Root cause

In `simulate_random_effects()` (`R/predictInterval_helpers.R`), the per-
level posterior variance is read as:

```r
reMatrix <- attr(rr[[j]], which = "postVar")
```

When the grouping factor `j` has a **single** correlated term block,
`postVar` is a 3-D array of shape `d × d × n_levels`. When the grouping
factor has **multiple** term blocks, lme4 returns `postVar` as a **list
of arrays**, each of shape `d_m × d_m × n_levels`, with one list element
per term block. The code downstream assumes an array and calls
`dimnames(reMatrix)[[3]] <- alllvl` (lines 270, 273) and indexes with
`reMatrix[, , k]` (line 284), which fails on a list.

Empirical verification:

| Formula                                        | `postVar` structure |
|------------------------------------------------|---------------------|
| `(Days + av  \| Subject)` (correlated)         | `array 3×3×18`      |
| `(Days + av \|\| Subject)` (uncorrelated)      | `list[3]` of `1×1×18` |
| `(Days\|Subject) + (0+av\|Subject)` (mixed)    | `list[2]` of `2×2×18` and `1×1×18` |

The combined dimension (sum of `d_m`) equals `ncol(reMeans)` in every
case, and the list ordering matches the column ordering of `reMeans`
and of `getME(merMod, "cnms")[[j]]`.

## Fix strategy

Normalize `postVar` to a 3-D block-diagonal array as soon as it is
read. This keeps the downstream code path — dimname assignment, level
filtering, `mvtnorm::rmvnorm()` — unchanged.

For each level `k`, the combined covariance matrix is block-diagonal
across term blocks (zero off-diagonals, because uncorrelated). For a
mixed case like `(Days|Subject) + (0+av|Subject)`, block 1 is a 2×2
full covariance matrix and block 2 is 1×1; the combined 3×3 matrix has
a 2×2 block and a 1×1 block on the diagonal, zeros elsewhere — exactly
the correct joint posterior covariance implied by the model
specification.

### Pseudocode

```r
reMatrix <- attr(rr[[j]], which = "postVar")
if (is.list(reMatrix)) {
  reMatrix <- combine_postvar_blocks(reMatrix, reMeans)
}
```

where `combine_postvar_blocks()` returns `array(dim = c(d, d, n))` with
block-diagonal structure per level, and carries the correct `dimnames`.

### Proposed helper (sketch)

```r
# Internal. Convert a per-term list of postVar arrays into a single
# block-diagonal array per level, matching the column structure of
# reMeans. Assumes list ordering matches colnames(reMeans), which lme4
# guarantees (verified via getME(, "cnms")).
combine_postvar_blocks <- function(pv_list, reMeans) {
  n_levels <- dim(pv_list[[1]])[3]
  d_total  <- ncol(reMeans)
  out <- array(0, dim = c(d_total, d_total, n_levels))
  pos <- 1L
  for (block in pv_list) {
    d_m <- dim(block)[1]
    idx <- seq.int(pos, length.out = d_m)
    for (lvl in seq_len(n_levels)) {
      out[idx, idx, lvl] <- block[, , lvl]
    }
    pos <- pos + d_m
  }
  dimnames(out) <- list(colnames(reMeans),
                        colnames(reMeans),
                        rownames(reMeans))
  out
}
```

Placement: private helper in `R/predictInterval_helpers.R`,
`@keywords internal`, no export.

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| lme4 may change list ordering in future versions | Add a runtime invariant: `sum(vapply(pv_list, function(b) dim(b)[1], 0L)) == ncol(reMeans)`. If it fails, stop with a clear message pointing to the lme4 version. |
| Different `d_m` per level (shouldn't happen) | The inner loop is per-level; shape is verified by `dim()` at the outer level. Not a practical concern. |
| Introducing non-zero off-diagonals the model didn't estimate | By construction we zero-initialize and only fill the diagonal blocks. `mvtnorm::rmvnorm()` will see the correct joint covariance. |
| Existing correlated-only models regress | `is.list(reMatrix)` guards the conversion; single-block (array) case is untouched. |
| GLMM `postVar` behavior differs | lme4 uses the same `postVar` attribute convention for GLMMs. Include one GLMM test case. |
| `mvtnorm::rmvnorm(method = "chol")` requires PD matrix; block-diag with a zero-variance block might be PSD-only | A zero-variance RE term is degenerate; lme4 would typically refit or warn before that. Still, add a test that asserts `sigma` is PD (all diagonal entries > 0) before the `rmvnorm` call, or fall back to `method = "svd"` if `chol` fails. Out of scope for this PR — document as a known edge in the NEWS entry. |

## Test coverage plan

All tests go in a new file `tests/testthat/test-predict-uncorrelated-re.R`
so Layer 1 (invariants) coverage grows without disturbing existing
suites. Each test uses `set.seed(11213)` per the repo convention.

1. **Regression: exact #118 reproducer**
   `lmer(Reaction ~ Days + (Days + av || Subject), data = sleepstudy)`
   Assert `predictInterval()` returns a 3-column `data.frame` with
   `nrow(newdata)` rows, finite values, and `upr >= fit >= lwr`.

2. **Double-bar with intercept**
   `lmer(Reaction ~ Days + (1 + Days || Subject), data = sleepstudy)`.
   Same structural assertions.

3. **Explicit split equivalent**
   `lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), ...)`.
   Should produce output equivalent in structure to test 2 (same
   dimensions, sensible intervals). Not bit-identical — different
   parameterization — but the marginal predictions should be within a
   tight tolerance.

4. **Mixed correlated + uncorrelated**
   `lmer(Reaction ~ Days + (Days|Subject) + (0 + av|Subject), ...)`.
   Exercises the 2×2 + 1×1 block-diagonal path.

5. **GLMM uncorrelated**
   Binomial with `(1 + x || g)` to confirm GLMM path is unaffected.

6. **No regression on correlated case**
   Existing snapshot/regression tests already cover this; no new test
   needed. Confirm the Layer 2 snapshot suite still passes unchanged.

7. **Helper unit test**
   Direct unit test on `combine_postvar_blocks()`: given a hand-built
   list of two arrays with known values, assert the combined array has
   the expected block structure and zero off-diagonals.

### Coverage guardrail

Before merge:
- `devtools::test()` must show ≥ previous pass count + 7 new tests
- `covr::file_coverage("R/predictInterval_helpers.R")` must not
  decrease

## Implementation checklist

- [ ] Add `combine_postvar_blocks()` helper in
      `R/predictInterval_helpers.R`
- [ ] Insert the `is.list(reMatrix)` branch at line ~258 to normalize
- [ ] Add runtime invariant check (sum of block dims == ncol(reMeans))
- [ ] Create `tests/testthat/test-predict-uncorrelated-re.R` with the
      seven tests above
- [ ] Run full suite locally (`devtools::test()`), confirm no
      regressions
- [ ] Run `R CMD check --as-cran` locally; must pass with 0 errors / 0
      warnings
- [ ] NEWS.md entry under "merTools 0.9.0 (unreleased) → Bug Fixes":
      one paragraph referencing #118
- [ ] Commit with message `fix: support uncorrelated random slopes in
      predictInterval (#118)`

## Out of scope

- Full audit of `predictInterval` semantics for unobserved levels in
  nested models — that is #124
- nlmer support — that is #121
- `REimpact` / `FEsim` interaction with uncorrelated terms — open a
  separate issue if a defect surfaces during testing
