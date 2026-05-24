# Combine a per-term list of postVar arrays into a block-diagonal array

When a grouping factor has multiple random-effect term blocks — from the
double-bar syntax (\`(x + y \|\| g)\`), explicit splits (\`(1\|g) +
(0+x\|g)\`), or mixed correlated + uncorrelated specs —
\`lme4::ranef(..., condVar = TRUE)\` returns \`postVar\` as a list of
arrays, one per term block. Each element has shape \`d_m × d_m ×
n_levels\`. This helper concatenates the blocks into a single \`d × d ×
n_levels\` array with block-diagonal structure per level (zeros
off-diagonal because uncorrelated). The block ordering and column
ordering match \`colnames(reMeans)\` / lme4's \`cnms\`.

## Usage

``` r
combine_postvar_blocks(pv_list, reMeans)
```

## Arguments

- pv_list:

  List of per-term postVar arrays.

- reMeans:

  The \`as.matrix(ranef(m)\[\[j\]\])\` matrix for the same grouping
  factor; used to size and name the combined array.

## Value

3-D array of shape \`ncol(reMeans) × ncol(reMeans) × nrow(reMeans)\`.
