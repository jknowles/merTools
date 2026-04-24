# Launch a shiny app to explore your merMod interactively

`shinyMer` launches a shiny app that allows you to interactively explore
an estimated merMod using functions from `merTools`.

## Usage

``` r
shinyMer(merMod, simData = NULL, pos = 1)
```

## Arguments

- merMod:

  An object of class "merMod".

- simData:

  A data.frame to make predictions from (optional). If NULL, then the
  user can only make predictions using the data in the frame slot of the
  merMod object.

- pos:

  The position of the environment to export function arguments to.
  Defaults to 1, the global environment, to allow shiny to run.

## Value

A shiny app
