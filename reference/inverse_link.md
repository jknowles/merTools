# Inverse link function for GLMM families

Internal helper that applies the inverse link function to convert linear
predictors to distribution parameters.

## Usage

``` r
inverse_link(eta, family, link)
```

## Arguments

- eta:

  Linear predictor (fixed + random effects).

- family:

  Character string specifying the distribution family.

- link:

  Character string specifying the link function.

## Value

Vector of distribution parameters (prob, lambda, or rate).
