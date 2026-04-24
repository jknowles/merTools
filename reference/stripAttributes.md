# Remove attributes from a data.frame

Strips attributes off of a data frame that come with a merMod
model.frame

## Usage

``` r
stripAttributes(data)
```

## Arguments

- data:

  a data.frame

## Value

a data frame with variable names cleaned to remove all attributes except
for names, row.names, and class
