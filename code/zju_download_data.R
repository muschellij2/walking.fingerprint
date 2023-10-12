library(purrr)
library(tidyverse)
# download data from https://www.ytzhang.net/datasets/zju-gaitacc/#citation
# put into data folder

read_text <- function(x) {
  read.table(x, sep = ",") %>% 
    t() %>% 
    as.data.frame() %>%
    mutate(
      ID = sub("/rec.*", "", sub(".*subj_", "", x)),
      session = sub("/subj.*", "", sub(".*gaitacc/", "", x)),
      record = sub(".txt.*", "", sub(".*rec_", "", x))
    )
}

files <-
  c(
    list.files(
      here::here("data/zju-gaitacc/session_0"),
      recursive = T,
      full.names = T
    ),
    list.files(
      here::here("data/zju-gaitacc/session_1"),
      recursive = T,
      full.names = T
    ),
    list.files(
      here::here("data/zju-gaitacc/session_2"),
      recursive = T,
      full.names = T
    )
  )


# remove names we don't need
vector <-
  c(grep("useful", files),
    grep("total", files),
    grep("cycles", files))
files <- files[-vector]

## read in all files to data frame
all_data <- files %>% map_dfr(read_text)

# get body locations and match IDs since session 0 individuals different from others
# t is the recording
all_data_zju_acc <-
  all_data %>%
  mutate(
    body_loc = sub(".*/", "", record),
    rec = sub("/.*", "", record),
    ID = ifelse(session == "session_0", as.numeric(ID), as.numeric(ID) + 22)
  ) %>% # to distinguish individuals
  group_by(ID, session, rec, body_loc) %>%
  mutate(t = row_number()) %>%
  ungroup()


# change your WD as necessary to get file list of useful files
files <-
  c(
    list.files(
      here::here("data/zju-gaitacc/session_0"),
      recursive = T,
      full.names = T
    ),
    list.files(
      here::here("data/zju-gaitacc/session_1"),
      recursive = T,
      full.names = T
    ),
    list.files(
      here::here("data/zju-gaitacc/session_2"),
      recursive = T,
      full.names = T
    )
  )

# remove names we don't need

vector <- c(grep("useful", files))
files <- files[vector]

# function to read in useful files
read_useful <- function(x) {
  read.table(x, sep = ",") %>% 
    rename("start" = V1, "end" = V2) %>%
    mutate(
      session = sub("/subj.*", "", sub(".*gaitacc/", "", x)),
      ID = sub("/rec.*", "", sub(".*subj_", "", x)),
      record = sub("/useful.*", "", sub(".*rec_", "", x))
    )
}

useful_key <- files %>% map_dfr(read_useful)
# now we have df with start and end of useful segments for each person, session, recording

# get IDS to match
useful_key <-
  useful_key %>%
  mutate(ID = ifelse(session == "session_0", as.numeric(ID), as.numeric(ID) + 22))


# left join with original data
joined <-
  all_data_zju_acc %>%
  left_join(.,
            useful_key,
            by = c(
              "session" = "session",
              "ID" = "ID",
              "rec" = "record"
            ))

# filter to just useful segment
filtered <-
  joined %>% 
  filter(t >= start & t <= end) %>% 
  dplyr::select(-c(start, end))

# change body location name
filtered <- filtered %>% mutate(
  loc = case_when(
    body_loc == "1" ~ "rwrist",
    body_loc == "2" ~ "larm",
    body_loc == "3" ~ "rpelvis",
    body_loc == "4" ~ "lthigh",
    body_loc == "5" ~ "rankle"
  ),
  signal = sqrt(V1 ^ 2 + V2 ^ 2 + V3 ^ 2)
) %>%
  dplyr::select(-c(body_loc, t)) %>%
  group_by(ID, session, loc) %>%
  mutate(s_allrec = row_number(),
         time_allrec = floor(s_allrec / 100) + 1) %>%
  ungroup() %>%
  mutate(second = ceiling(s_allrec / 100)) %>%
  rename(time = s_allrec) %>% 
  dplyr::select(ID, session, rec, loc, signal, second, time)


filtered_wide <-
  filtered %>% pivot_wider(names_from = loc,
                           values_from = signal,
                           names_prefix = "signal_")

write.csv(filtered_wide, here::here("data/df_all_zju.csv"))

rm(list = ls())
