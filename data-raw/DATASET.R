## code to prepare `DATASET` dataset goes here
require(dplyr)
experiment1 <- 
  read.table('./data-raw/Burr_2017_PDL1.txt') %>%
  head()
usethis::use_data(experiment1)
