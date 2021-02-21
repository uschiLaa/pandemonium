library(tidyverse)

b_anomaly <- list()

b_anomaly$pred <- read_csv("data-raw/pred.csv", col_names = F)
b_anomaly$wc <- read_csv("data-raw/wc_grid.csv", col_names = c("C9", "C10"))
b_anomaly$covInv <- read_csv("data-raw/cov_inv.csv", col_names = F) %>%
  as.matrix()
b_anomaly$exp <- read_csv("data-raw/exp.csv", col_names = c("value"))

usethis::use_data(b_anomaly)
