# Collapse a dataframe to a single average row

Take an entire dataframe and summarize it in one row by using the mean
and mode.

## Usage

``` r
collapseFrame(data)
```

## Arguments

- data:

  a data.frame

## Value

a data frame with a single row

## Details

Each character and factor variable in the data.frame is assigned to the
modal category and each numeric variable is collapsed to the mean.
Currently if mode is a tie, returns a "."
