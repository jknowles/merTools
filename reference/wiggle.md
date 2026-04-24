# Assign an observation to different values

Creates a new data.frame with copies of the original observation, each
assigned to a different user-specified value of a variable. Allows the
user to look at the effect on predicted values of changing either a
single variable or multiple variables.

## Usage

``` r
wiggle(data, varlist, valueslist)
```

## Arguments

- data:

  a data frame with one or more observations to be reassigned

- varlist:

  a character vector specifying the name(s) of the variable to adjust

- valueslist:

  a list of vectors with the values to assign to var

## Value

a `data.frame` with each row assigned to the one of the new variable
combinations. All variable combinations are returned, eg wiggling two
variables with 3 and 4 variables respectively will return a new dataset
with `3 * 4 = 12` observations.

## Details

If the variable specified is a factor, then wiggle will return it as a
character.

## Examples

``` r
data(iris)
wiggle(iris[3,], varlist = "Sepal.Width", valueslist = list(c(1, 2, 3, 5)))
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 3           4.7           1          1.3         0.2  setosa
#> 31          4.7           2          1.3         0.2  setosa
#> 32          4.7           3          1.3         0.2  setosa
#> 33          4.7           5          1.3         0.2  setosa
wiggle(iris[3:5,], "Sepal.Width", valueslist = list(c(1, 2, 3, 5)))
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 3           4.7           1          1.3         0.2  setosa
#> 4           4.6           1          1.5         0.2  setosa
#> 5           5.0           1          1.4         0.2  setosa
#> 31          4.7           2          1.3         0.2  setosa
#> 41          4.6           2          1.5         0.2  setosa
#> 51          5.0           2          1.4         0.2  setosa
#> 32          4.7           3          1.3         0.2  setosa
#> 42          4.6           3          1.5         0.2  setosa
#> 52          5.0           3          1.4         0.2  setosa
#> 33          4.7           5          1.3         0.2  setosa
#> 43          4.6           5          1.5         0.2  setosa
#> 53          5.0           5          1.4         0.2  setosa
wiggle(iris[3,], c("Sepal.Width", "Petal.Length"), list(c(1,2,3,5), c(3,5,6)))
#>                 Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> Sepal.Width.3            4.7           1            3         0.2  setosa
#> Sepal.Width.31           4.7           2            3         0.2  setosa
#> Sepal.Width.32           4.7           3            3         0.2  setosa
#> Sepal.Width.33           4.7           5            3         0.2  setosa
#> Sepal.Width.34           4.7           1            5         0.2  setosa
#> Sepal.Width.311          4.7           2            5         0.2  setosa
#> Sepal.Width.321          4.7           3            5         0.2  setosa
#> Sepal.Width.331          4.7           5            5         0.2  setosa
#> Sepal.Width.35           4.7           1            6         0.2  setosa
#> Sepal.Width.312          4.7           2            6         0.2  setosa
#> Sepal.Width.322          4.7           3            6         0.2  setosa
#> Sepal.Width.332          4.7           5            6         0.2  setosa
wiggle(iris[3:5,], c("Sepal.Width", "Petal.Length"), list(c(1,2,3,5), c(3,5,6)))
#>                 Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> Sepal.Width.3            4.7           1            3         0.2  setosa
#> Sepal.Width.4            4.6           1            3         0.2  setosa
#> Sepal.Width.5            5.0           1            3         0.2  setosa
#> Sepal.Width.31           4.7           2            3         0.2  setosa
#> Sepal.Width.41           4.6           2            3         0.2  setosa
#> Sepal.Width.51           5.0           2            3         0.2  setosa
#> Sepal.Width.32           4.7           3            3         0.2  setosa
#> Sepal.Width.42           4.6           3            3         0.2  setosa
#> Sepal.Width.52           5.0           3            3         0.2  setosa
#> Sepal.Width.33           4.7           5            3         0.2  setosa
#> Sepal.Width.43           4.6           5            3         0.2  setosa
#> Sepal.Width.53           5.0           5            3         0.2  setosa
#> Sepal.Width.34           4.7           1            5         0.2  setosa
#> Sepal.Width.44           4.6           1            5         0.2  setosa
#> Sepal.Width.54           5.0           1            5         0.2  setosa
#> Sepal.Width.311          4.7           2            5         0.2  setosa
#> Sepal.Width.411          4.6           2            5         0.2  setosa
#> Sepal.Width.511          5.0           2            5         0.2  setosa
#> Sepal.Width.321          4.7           3            5         0.2  setosa
#> Sepal.Width.421          4.6           3            5         0.2  setosa
#> Sepal.Width.521          5.0           3            5         0.2  setosa
#> Sepal.Width.331          4.7           5            5         0.2  setosa
#> Sepal.Width.431          4.6           5            5         0.2  setosa
#> Sepal.Width.531          5.0           5            5         0.2  setosa
#> Sepal.Width.35           4.7           1            6         0.2  setosa
#> Sepal.Width.45           4.6           1            6         0.2  setosa
#> Sepal.Width.55           5.0           1            6         0.2  setosa
#> Sepal.Width.312          4.7           2            6         0.2  setosa
#> Sepal.Width.412          4.6           2            6         0.2  setosa
#> Sepal.Width.512          5.0           2            6         0.2  setosa
#> Sepal.Width.322          4.7           3            6         0.2  setosa
#> Sepal.Width.422          4.6           3            6         0.2  setosa
#> Sepal.Width.522          5.0           3            6         0.2  setosa
#> Sepal.Width.332          4.7           5            6         0.2  setosa
#> Sepal.Width.432          4.6           5            6         0.2  setosa
#> Sepal.Width.532          5.0           5            6         0.2  setosa
```
