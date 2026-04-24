# Create a factor with unobserved levels

Create a factor variable and include unobserved levels for compatibility
with model prediction functions

## Usage

``` r
superFactor(x, fullLev)
```

## Arguments

- x:

  a vector to be converted to a factor

- fullLev:

  a vector of factor levels to be assigned to x

## Value

a factor variable with all observed levels of x and all levels of x in
fullLev

## Examples

``` r
# \donttest{
regularFactor <- c("A", "B", "C")
regularFactor <- factor(regularFactor)
levels(regularFactor)
#> [1] "A" "B" "C"
# Now make it super
newLevs <- c("D", "E", "F")
regularFactor <- superFactor(regularFactor, fullLev = newLevs)
levels(regularFactor) # now super
#> [1] "D" "E" "F" "A" "B" "C"
# }
```
