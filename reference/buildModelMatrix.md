# Build model matrix

a function to create a model matrix with all predictor terms in both the
group level and fixed effect level

## Usage

``` r
buildModelMatrix(model, newdata, which = "full")
```

## Source

Taken from predict.merMod in lme4

## Arguments

- model:

  a merMod object from lme4

- newdata:

  a data frame to construct the matrix from

- which:

  a character which matrix to return,default is full matrix with fixed
  and random terms, other options are "fixed" and "random"
