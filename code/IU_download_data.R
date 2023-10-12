library(here)
library(purrr)
library(tidyverse)

# create data folder
if (!dir.exists("data")) {
  dir.create(here::here("data"))
}

# get files
system(
  "wget -r -N -c -np https://physionet.org/files/accelerometry-walk-climb-drive/1.0.0/"
)

# get names of files
files <-
  list.files(
    here(
      "physionet.org/files/accelerometry-walk-climb-drive/1.0.0/raw_accelerometry_data"
    ),
    full.names = T
  )


# function to read files and create ID column
read_files <- function(x) {
  readr::read_csv(x, show_col_types = F) %>% mutate(ID = sub(".csv.*", "", sub(".*id", "", x)))
}

# read in all files and row bind into df
df_all <-
  files %>%
  map_dfr(read_files) %>%
  dplyr::select(-`<html>`)

df_all <- files %>% map_dfr(read_files)


df_all_filtered <-
  df_all  %>%
  filter(activity == 1) %>% # filtered to just walking
  mutate(time_g = time_s - min(time_s), # starting time from 0
         time_g = floor(time_g)) %>%
  group_by(ID) %>%
  mutate(time_g = time_g - min(time_g)) %>%
  group_by(ID, time_g) %>% # group by ID and second, get number of observations per second
  mutate(n_g = n(),
         time_s_g = 1:n()) %>%
  ungroup() %>%
  filter(n_g == 100) %>% # keeping just seconds with whole second of observations
  dplyr::select(-n_g) %>%
  mutate(
    signal_lw = sqrt(lw_x ^ 2 + lw_y ^ 2 + lw_z ^ 2),
    signal_lh = sqrt(lh_x ^ 2 + lh_y ^ 2 + lh_z ^ 2),
    signal_la = sqrt(la_x ^ 2 + la_y ^ 2 + la_z ^ 2),
    signal_ra = sqrt(ra_x ^ 2 + ra_y ^ 2 + ra_z ^ 2),
  ) %>% mutate(ID2 = as.numeric(as.factor(ID))) # making ID numeric

df_all_final <-
  df_all_filtered %>%
  group_by(ID2) %>%
  mutate(rownum = row_number(),
         time_s_2 = rownum / 100,
         J = ceiling(time_s_2)) %>%
  ungroup()


## get columns in unified format
df_all_final <-
  df_all_final %>%
  dplyr::select(signal_lw, signal_lh, signal_la, signal_ra, ID2, time_s_2, J) %>%
  rename(time = time_s_2,
         second = J)

write.csv(df_all_final, here::here("data/df_all_IU.csv"))

rm(list = ls())


