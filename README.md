
# ageage

<!-- badges: start -->
<!-- badges: end -->

example R/Rcpp package

## Installation

``` r
# install.packages("devtools")
devtools::install_github("ben-williams/ageage")
```

## Example


``` r
# load ----
library(ageage)

# analysis ----
# setup folders
setup(year = 2022)

# age error for GOA POP
# if your admb lives at "c:/admb" can leave NULL, othwerwise point it to the correct place
age_error(year = 2022, admb_home = NULL, rec_age = 2, plus_age = 25)

```

