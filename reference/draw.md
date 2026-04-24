# Draw a single observation out of an object matching some criteria

Draw is used to select a single observation out of an R object.
Additional parameters allow the user to control how that observation is
chosen in order to manipulate that observation later. This is a generic
function with methods for a number of objects.

## Usage

``` r
draw(object, type = c("random", "average"), varList = NULL, seed = NULL, ...)

# S3 method for class 'merMod'
draw(object, type = c("random", "average"), varList = NULL, seed = NULL, ...)
```

## Arguments

- object:

  the object to draw from

- type:

  what kind of draw to make. Options include random or average

- varList:

  a list specifying filters to subset the data by when making the draw

- seed:

  numeric, optional argument to set seed for simulations, ignored if
  type="average"

- ...:

  additional arguments required by certain methods

## Value

a data.frame with a single row representing the desired observation

## Details

In cases of tie, ".", may be substituted for factors.

## Examples

``` r
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
# Random case
draw(fm1, type = "random")
#>    Reaction Days Subject
#> 27 235.4208    6     310
# Average
draw(fm1, type = "average")
#> Number of observations < 20, random effect quantiles may not be well-defined.
#>   Reaction Days Subject
#> 1 298.5079  4.5     308
# Subset
draw(fm1, type = "average", varList = list("Subject" = "308"))
#> Warning: Subset has less than 20 rows, averages may be problematic.
#>   Reaction Days Subject
#> 1 342.1338  4.5     308
```
