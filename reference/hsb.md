# A subset of data from the 1982 High School and Beyond survey used as examples for HLM software

A key example dataset used for examples in the HLM software manual.
Included here for use in replicating HLM analyses in R.

## Usage

``` r
hsb
```

## Format

A data frame with 7,185 observations on the following 8 variables.

- `schid`:

  a numeric vector, 160 unique values

- `mathach`:

  a numeric vector for the performance on a standardized math assessment

- `female`:

  a numeric vector coded 0 for male and 1 for female

- `ses`:

  a numeric measure of student socio-economic status

- `minority`:

  a numeric vector coded 0 for white and 1 for non-white students

- `schtype`:

  a numeric vector coded 0 for public and 1 for private schools

- `meanses`:

  a numeric, the average SES for each school in the data set

- `size`:

  a numeric for the number of students in the school

## Source

Data made available by UCLA Institute for Digital Research and Education
(IDRE) online:
<https://stats.oarc.ucla.edu/other/hlm/hlm-mlm/introduction-to-multilevel-modeling-using-hlm/>

## Details

The data file used for this presentation is a subsample from the 1982
High School and Beyond Survey and is used extensively in Hierarchical
Linear Models by Raudenbush and Bryk. It consists of 7,185 students
nested in 160 schools.

## References

Stephen W. Raudenbush and Anthony S. Bryk (2002). Hierarchical Linear
Models: Applications and Data Analysis Methods (2nd ed.). SAGE.

## Examples

``` r
data(hsb)
head(hsb)
#>   schid minority female    ses mathach size schtype meanses
#> 1  1224        0      1 -1.528   5.876  842       0  -0.428
#> 2  1224        0      1 -0.588  19.708  842       0  -0.428
#> 3  1224        0      0 -0.528  20.349  842       0  -0.428
#> 4  1224        0      0 -0.668   8.781  842       0  -0.428
#> 5  1224        0      0 -0.158  17.898  842       0  -0.428
#> 6  1224        0      0  0.022   4.583  842       0  -0.428
```
