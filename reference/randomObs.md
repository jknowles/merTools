# Select a random observation from model data

Select a random observation from the model frame of a merMod

## Usage

``` r
randomObs(merMod, varList, seed = NULL)
```

## Arguments

- merMod:

  an object of class merMod

- varList:

  optional, a named list of conditions to subset the data on

- seed:

  numeric, optional argument to set seed for simulations

## Value

a data frame with a single row for a random observation, but with full
factor levels. See details for more.

## Details

Each factor variable in the data frame has all factor levels from the
full model.frame stored so that the new data is compatible with
predict.merMod
