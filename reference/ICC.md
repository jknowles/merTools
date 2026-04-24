# Calculate the intraclass correlation using mixed effect models

Calculate the intraclass correlation using mixed effect models

## Usage

``` r
ICC(outcome, group, data, subset = NULL)
```

## Arguments

- outcome:

  a character representing the variable of the outcome

- group:

  a character representing the name of the grouping term

- data:

  a data.frame

- subset:

  an optional subset

## Value

a numeric for the intraclass correlation

## Examples

``` r
data(sleepstudy)
ICC(outcome = "Reaction", group = "Subject", data = sleepstudy)
#> [1] 0.3948896
```
