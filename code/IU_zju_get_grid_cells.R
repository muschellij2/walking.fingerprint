library(readr)
library(purrr)
library(tidyverse)

# this file gets the number of points in each grid cell for both datasets 

df_all_zju <- read_csv(here::here("data/df_all_zju.csv"),
                       col_types = cols(...1 = col_skip()))

df_all_IU <- read_csv(here::here("data/df_all_IU.csv"),
                      col_types = cols(...1 = col_skip())) %>%
  rename(
    ID = ID2,
    signal_lwrist = signal_lw,
    signal_lhip = signal_lh,
    signal_lankle = signal_la,
    signal_rankle = signal_ra
  )



# function to get grid cell data from one subject 
get_grid_data  <- function(subject, time_lags, gcell_size, location, data){
  # function to get grid cell data from one subject
  # select ID, signal, time (seconds), filter to subject
  df <-
    data %>%
    dplyr::select(ID, paste("signal_", location, sep = ""), second) %>%
    rename(signal = paste("signal_", location, sep = "")) %>%
    filter(ID == subject)
  
  # max_signal <- round(max(df$signal), 0)
  max_signal <-
    3  # we set max signal to 3 based on EDA, but could take actual max signal
  n_seconds <-
    max(df$second) # number of total seconds for the subject
  seconds <-
    rep(seq(1, n_seconds, 1), each = length(time_lags)) # vector of seconds and lags so that we can iterate over both
  lags <- rep(time_lags, n_seconds) # vector of lags
  # function to get the grid data for one second, one lag
  get_grid_data_lagsec <- function(s, lag) {
    # filter to one second
    df %>% filter(second == s) %>%
      dplyr::select(signal) %>%
      mutate(lag_signal = dplyr::lag(signal, n = lag)) %>%   # for each second, calculate signal and lagged signal
      mutate(
        cut_sig = cut(
          signal,
          breaks = seq(0, max_signal, by = gcell_size),
          include.lowest = T
        ),
        cut_lagsig = cut(
          lag_signal,
          breaks = seq(0, max_signal, by = gcell_size),
          include.lowest = T
        )
      ) %>%
      drop_na() %>% # count # points in each "grid cell"
      count(cut_sig, cut_lagsig, .drop = FALSE) %>%
      mutate(
        lag_hz = lag,
        ID = subject,
        second = s,
        cell = paste(cut_sig, cut_lagsig, lag_hz)
      ) %>%
      dplyr::select(n, ID, second, cell)
  }
  # apply above function over all seconds and lags, return df
  map2_dfr(.x = seconds, .y = lags,
           .f = get_grid_data_lagsec) %>% pivot_wider(
             id_cols = c(ID, second),
             names_from = cell,
             values_from = n
           )
  # apply across all seconds and lags
}

subs_IU <- unique(df_all_IU$ID)

# get left wrist data for IU 
grid_data_lw_IU <- map_dfr(.x = subs_IU,
                           .f = get_grid_data,
                           time_lags = c(15, 30, 45),
                           gcell_size = 0.25,
                           location = "lwrist",
                           data = df_all_IU,
                           .progress = T)
saveRDS(grid_data_lw_IU, file = here::here("data/grid_data_lw_IU.rds"))

subs_zju <- unique(df_all_zju$ID[df_all_zju$session!="session_0"]) 

grid_data_rw_zju_s1 <- map_dfr(.x = subs_zju,
                               .f = get_grid_data,
                               time_lags = c(15, 30, 45),
                               gcell_size = 0.25,
                               location = "rwrist",
                               data = df_all_zju %>% filter(session=="session_1"),
                               .progress = T)

grid_data_rw_zju_s2 <- map_dfr(.x = subs_zju,
                               .f = get_grid_data,
                               time_lags = c(15, 30, 45),
                               gcell_size = 0.25,
                               location = "rwrist",
                               data = df_all_zju %>% filter(session=="session_2"),
                               .progress = T)

saveRDS(grid_data_rw_zju_s1, file = here::here("data/grid_data_rw_zju_s1.rds"))
saveRDS(grid_data_rw_zju_s2, file = here::here("data/grid_data_rw_zju_s2.rds"))


rm(list=ls())
