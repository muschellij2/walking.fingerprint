## code to prepare `example_walking` dataset goes here
library(dplyr)
library(adeptdata)
example_walking = adeptdata::acc_walking_IU %>%
  filter(subj_id == subj_id[1], loc_id == "left_wrist") %>%
  mutate(vm = sqrt(x^2 + y^2 + z^2)) %>%
  select(time = time_s, vm)
example_walking = dplyr::as_tibble(example_walking)
usethis::use_data(example_walking, overwrite = TRUE)
