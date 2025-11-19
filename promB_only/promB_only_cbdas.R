setwd("C:/Users/scej4/OneDrive/Desktop/diplomska_koda/promB_only")
library(tidyverse)

plus_strand <- function(x) {
  read.table(file = x, sep = "\t") |>
    filter(V6 == '+')
}

minus_strand <- function(x) {
  read.table(file = x, sep = "\t") |>
    filter(V6 == '-')
}
