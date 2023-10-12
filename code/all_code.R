# IU_download_data 

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

# zju_download_data 
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

# IU_ZJU_get_grid_cells

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

# IU_fit_logistic
library(tidyverse)
library(purrr)
library(tidymodels)
library(magrittr)
tidymodels_prefer()
options(dplyr.summarise.inform = FALSE)

## read in data
grid_data_lw_IU <- readRDS(here::here("data/grid_data_lw_IU.rds"))

# split into training and testing
# 75% train, 25% test, equal proportions for ea person
data_split <- split(grid_data_lw_IU, f = grid_data_lw_IU$ID)
# function to sample percentage
samp <- function(pct, n, ind) {
  set.seed(ind)
  sample(n, floor(pct * n), replace = F)
}
ids <- unique(grid_data_lw_IU$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

getrows <- function(data, rows) {
  data[rows, ]
}
getrows_test <- function(data, rows) {
  data[-rows, ]
}

data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)

# function to fit ovr logistic regression models
fit_model <- function(subject, train, test) {
  train$class <- ifelse(train$ID == subject, 1, 0)
  test$class <- ifelse(test$ID == subject, 1, 0)
  tmp <- train %>% dplyr::select(-c(ID, second))
  tmp_test <- test %>% dplyr::select(-c(ID, second))
  mod <-
    glm(class ~ ., data = tmp, family = binomial(link = "logit"))
  pred <- predict.glm(mod, newdata = tmp_test, type = "response")
  return(pred)
}

# first we want to remove columns with near zero variance
nzv_trans <-
  recipe(ID ~ ., data = data_train) %>%
  step_nzv(all_predictors())

nzv_estimates <- prep(nzv_trans)

nzv <- colnames(juice(nzv_estimates))
dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv))
dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv))

## now fit models, get predictions
t1 <- Sys.time()
all_predictions_IU <-
  map_dfc(
    .x = ids,
    .f = fit_model,
    train = dat_nzv,
    test = dat_nzv_test,
    .progress = T
  ) %>%
  janitor::clean_names()

Sys.time() - t1
# Time difference of 44.60639 secs


# column j is predicted probability that data in that row belong to subject j
# normalize probabilities
row_sums <- rowSums(all_predictions_IU)

# normalize and add "true subject column"
all_predictions_IU %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x32, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = dat_nzv_test$ID)



if (!dir.exists("predictions")) {
  dir.create(here::here("predictions"))
}
saveRDS(all_predictions_IU,
        here::here("predictions/IU_logistic_predictions.rds"))

# print some results summaries

get_summarized_predictions <- function(predictions, long = FALSE) {
  if (long == T) {
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      mutate(correct = ifelse(true_subject == model, 1, 0))
  }
  else{
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      group_by(true_subject) %>%
      summarize(
        maxprob = first(max(mean_pred)),
        predicted_sub = first(model[mean_pred == maxprob]),
        probsubj = first(mean_pred[true_subject == model])
      ) %>%
      mutate(correct = ifelse(as.numeric(predicted_sub) == true_subject, 1, 0))
  }
}

get_summarized_predictions(all_predictions_IU) %>% print(n = Inf)

# fn takes number of seconds to average over as input, outputs classification stats
get_prediction_stats <- function(predictions, seconds) {
  predictions %>%
    group_by(true_subject) %>%
    mutate(sec = floor(row_number() / seconds)) %>%
    pivot_longer(cols = -c("true_subject", "sec")) %>%
    mutate(model = as.numeric(sub(".*x", "", name))) %>%
    rename(pred = value) %>%
    ungroup() %>%
    group_by(true_subject, model, sec) %>%
    summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
    group_by(true_subject, sec) %>%
    summarize(maxprob = max(mean_pred),
              predicted_subject = model[mean_pred == maxprob]) %>% ungroup() %>%
    dplyr::select(c(true_subject, predicted_subject)) %>%
    mutate(across(1:2, as.factor)) %>%
    yardstick::conf_mat(., truth = true_subject, estimate = predicted_subject) %>%
    summary() %>%
    mutate(s = seconds)
}

results_summarized <-
  map_dfr(
    .x = c(1, seq(10, 100, 10)),
    .f = get_prediction_stats,
    predictions = all_predictions_IU,
    .progress = T
  ) %>%
  filter(.metric != "detection_prevalence")
supp.labs <-
  c(
    "Accuracy",
    "Kappa",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NVP",
    "MCC",
    "J Index",
    "Balanced Accuracy",
    "Precision",
    "Recall",
    "F1 Score"
  )
names(supp.labs) <- unique(results_summarized$.metric)

results_summarized %>%
  ggplot(aes(x = s, y = .estimate, col = .metric)) +
  geom_point() +
  geom_line() +
  theme_light() +
  facet_wrap(. ~ .metric,
             scales = "free_y",
             labeller = labeller(.metric = supp.labs)) +
  labs(x = "Number of Seconds", y = "Estimate") +
  scale_x_continuous(breaks = c(1, seq(10, 100, 10))) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "none")

rm(list = ls())

# zju_fit_logistic
library(tidyverse)
library(purrr)
library(tidymodels)
library(magrittr)
tidymodels_prefer()
options(dplyr.summarise.inform = FALSE)

## read in data
grid_data_rw_zju_s1 <-
  readRDS(here::here("data/grid_data_rw_zju_s1.rds")) %>%
  mutate(ID = ID - 22)
grid_data_rw_zju_s2 <-
  readRDS(here::here("data/grid_data_rw_zju_s2.rds")) %>%
  mutate(ID = ID - 22)

# first: use just session 1 to predict on session 1
# split into training and testing
# 75% train, 25% test, equal proportions for ea person
data_split <- split(grid_data_rw_zju_s1, f = grid_data_rw_zju_s1$ID)
# function to sample percentage
samp <- function(pct, n, ind) {
  set.seed(ind)
  sample(n, floor(pct * n), replace = F)
}
ids <- unique(grid_data_rw_zju_s1$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

getrows <- function(data, rows) {
  data[rows,]
}
getrows_test <- function(data, rows) {
  data[-rows,]
}

data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)

# function to fit ovr logistic regression models
fit_model <- function(subject, train, test) {
  train$class <- ifelse(train$ID == subject, 1, 0)
  test$class <- ifelse(test$ID == subject, 1, 0)
  tmp <- train %>% dplyr::select(-c(ID, second))
  tmp_test <- test %>% dplyr::select(-c(ID, second))
  mod <-
    glm(class ~ ., data = tmp, family = binomial(link = "logit"))
  pred <- predict.glm(mod, newdata = tmp_test, type = "response")
  return(pred)
}

# first we want to remove columns with near zero variance
nzv_trans <-
  recipe(ID ~ ., data = data_train) %>%
  step_nzv(all_predictors())

nzv_estimates <- prep(nzv_trans)

nzv <- colnames(juice(nzv_estimates))
dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv))
dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv))

## now fit models, get predictions
t1 <- Sys.time()
all_predictions_zju_s1 <-
  map_dfc(
    .x = ids,
    .f = fit_model,
    train = dat_nzv,
    test = dat_nzv_test,
    .progress = T
  ) %>%
  janitor::clean_names()
Sys.time() - t1

# Time difference of 2.041705 mins

row_sums <- rowSums(all_predictions_zju_s1)

# normalize and add "true subject column"
all_predictions_zju_s1 %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = dat_nzv_test$ID)



if (!dir.exists("predictions")) {
  dir.create(here::here("predictions"))
}
saveRDS(
  all_predictions_zju_s1,
  here::here("predictions/zjus1_logistic_predictions.rds")
)

get_summarized_predictions <- function(predictions, long = FALSE) {
  if (long == T) {
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      mutate(correct = ifelse(true_subject == model, 1, 0))
  }
  else{
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      group_by(true_subject) %>%
      summarize(
        maxprob = first(max(mean_pred)),
        predicted_sub = first(model[mean_pred == maxprob]),
        probsubj = first(mean_pred[true_subject == model])
      ) %>%
      mutate(correct = ifelse(as.numeric(predicted_sub) == true_subject, 1, 0))
  }
}

get_summarized_predictions(all_predictions_zju_s1) %>% print(n = Inf)

get_prediction_stats <- function(predictions, seconds) {
  predictions %>%
    group_by(true_subject) %>%
    mutate(sec = floor(row_number() / seconds)) %>%
    pivot_longer(cols = -c("true_subject", "sec")) %>%
    mutate(model = as.numeric(sub(".*x", "", name))) %>%
    rename(pred = value) %>%
    ungroup() %>%
    group_by(true_subject, model, sec) %>%
    summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
    group_by(true_subject, sec) %>%
    summarize(maxprob = max(mean_pred),
              predicted_subject = model[mean_pred == maxprob]) %>% ungroup() %>%
    dplyr::select(c(true_subject, predicted_subject)) %>%
    mutate(across(1:2, factor, levels = as.factor(seq(1, 153, 1)))) %>%
    yardstick::conf_mat(., truth = true_subject, estimate = predicted_subject) %>%
    summary() %>%
    mutate(s = seconds)
}

results_summarized <-
  map_dfr(
    .x = seq(1, 21, 2),
    .f = get_prediction_stats,
    predictions = all_predictions_zju_s1,
    .progress = T
  ) %>%
  filter(.metric != "detection_prevalence")

supp.labs <-
  c(
    "Accuracy",
    "Kappa",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NVP",
    "MCC",
    "J Index",
    "Balanced Accuracy",
    "Precision",
    "Recall",
    "F1 Score"
  )
names(supp.labs) <- unique(results_summarized$.metric)

results_summarized %>%
  ggplot(aes(x = s, y = .estimate, col = .metric)) +
  geom_point() +
  geom_line() +
  theme_light() +
  facet_wrap(. ~ .metric,
             scales = "free_y",
             labeller = labeller(.metric = supp.labs)) +
  labs(x = "Number of Seconds", y = "Estimate") +
  scale_x_continuous(breaks = seq(1, 21, 2)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "none")


## train on session 2, predict on session 2
data_split <- split(grid_data_rw_zju_s2, f = grid_data_rw_zju_s2$ID)

ids <- unique(grid_data_rw_zju_s2$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)

# function to fit ovr logistic regression models

# first we want to remove columns with near zero variance
nzv_trans <-
  recipe(ID ~ ., data = data_train) %>%
  step_nzv(all_predictors())

nzv_estimates <- prep(nzv_trans)

nzv <- colnames(juice(nzv_estimates))
dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv))
dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv))

## now fit models, get predictions
t1 <- Sys.time()
all_predictions_zju_s2 <-
  map_dfc(
    .x = ids,
    .f = fit_model,
    train = dat_nzv,
    test = dat_nzv_test,
    .progress = T
  ) %>%
  janitor::clean_names()
Sys.time() - t1

# Time difference of 2.299532 mins

row_sums <- rowSums(all_predictions_zju_s2)

# normalize and add "true subject column"
all_predictions_zju_s2 %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = dat_nzv_test$ID)


saveRDS(
  all_predictions_zju_s2,
  here::here("predictions/zjus2_logistic_predictions.rds")
)

get_summarized_predictions(all_predictions_zju_s2) %>% print(n = Inf)

results_summarized <-
  map_dfr(
    .x = seq(1, 21, 2),
    .f = get_prediction_stats,
    predictions = all_predictions_zju_s2,
    .progress = T
  ) %>%
  filter(.metric != "detection_prevalence")

supp.labs <-
  c(
    "Accuracy",
    "Kappa",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NVP",
    "MCC",
    "J Index",
    "Balanced Accuracy",
    "Precision",
    "Recall",
    "F1 Score"
  )
names(supp.labs) <- unique(results_summarized$.metric)

results_summarized %>%
  ggplot(aes(x = s, y = .estimate, col = .metric)) +
  geom_point() +
  geom_line() +
  theme_light() +
  facet_wrap(. ~ .metric,
             scales = "free_y",
             labeller = labeller(.metric = supp.labs)) +
  labs(x = "Number of Seconds", y = "Estimate") +
  scale_x_continuous(breaks = seq(1, 21, 2)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "none")



# all data combined
grid_data_all <-
  bind_rows(grid_data_rw_zju_s1, grid_data_rw_zju_s2) %>%
  group_by(ID) %>%
  mutate(second = row_number())

data_split <- split(grid_data_all, f = grid_data_all$ID)

ids <- unique(grid_data_all$ID)
# number of rows for each individual
data_split <- split(grid_data_all, f = grid_data_all$ID)


ids <- unique(grid_data_all$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)

# function to fit ovr logistic regression models

# first we want to remove columns with near zero variance
nzv_trans <-
  recipe(ID ~ ., data = data_train) %>%
  step_nzv(all_predictors())

nzv_estimates <- prep(nzv_trans)

nzv <- colnames(juice(nzv_estimates))
dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv))
dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv))

## now fit models, get predictions
t1 <- Sys.time()
all_predictions_zju_all <-
  map_dfc(
    .x = ids,
    .f = fit_model,
    train = dat_nzv,
    test = dat_nzv_test,
    .progress = T
  ) %>%
  janitor::clean_names()
Sys.time() - t1

# Time difference of 4.015318 mins

row_sums <- rowSums(all_predictions_zju_all)

# normalize and add "true subject column"
all_predictions_zju_all %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = dat_nzv_test$ID)



saveRDS(
  all_predictions_zju_all,
  here::here("predictions/zjuall_logistic_predictions.rds")
)

get_summarized_predictions(all_predictions_zju_all) %>% print(n = Inf)

results_summarized <-
  map_dfr(
    .x = seq(1, 21, 2),
    .f = get_prediction_stats,
    predictions = all_predictions_zju_all,
    .progress = T
  ) %>%
  filter(.metric != "detection_prevalence")

supp.labs <-
  c(
    "Accuracy",
    "Kappa",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NVP",
    "MCC",
    "J Index",
    "Balanced Accuracy",
    "Precision",
    "Recall",
    "F1 Score"
  )
names(supp.labs) <- unique(results_summarized$.metric)

results_summarized %>%
  ggplot(aes(x = s, y = .estimate, col = .metric)) +
  geom_point() +
  geom_line() +
  theme_light() +
  facet_wrap(. ~ .metric,
             scales = "free_y",
             labeller = labeller(.metric = supp.labs)) +
  labs(x = "Number of Seconds", y = "Estimate") +
  scale_x_continuous(breaks = seq(1, 21, 2)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "none")

# finally: train on 1, predict on 2
data_train <- grid_data_rw_zju_s1
data_test <- grid_data_rw_zju_s2


# first we want to remove columns with near zero variance
nzv_trans <-
  recipe(ID ~ ., data = data_train) %>%
  step_nzv(all_predictors())

nzv_estimates <- prep(nzv_trans)

nzv <- colnames(juice(nzv_estimates))
dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv))
dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv))

## now fit models, get predictions
t1 <- Sys.time()
all_predictions_zju_s1s2 <-
  map_dfc(
    .x = ids,
    .f = fit_model,
    train = dat_nzv,
    test = dat_nzv_test,
    .progress = T
  ) %>%
  janitor::clean_names()

Sys.time() - t1

# Time difference of 2.71424 mins

row_sums <- rowSums(all_predictions_zju_s1s2)

# normalize and add "true subject column"
all_predictions_zju_s1s2 %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = dat_nzv_test$ID)


saveRDS(
  all_predictions_zju_s1s2,
  here::here("predictions/zjus1s2_logistic_predictions.rds")
)

get_summarized_predictions(all_predictions_zju_s1s2) %>% print(n = Inf)

results_summarized <-
  map_dfr(
    .x = seq(1, 21, 2),
    .f = get_prediction_stats,
    predictions = all_predictions_zju_s1s2,
    .progress = T
  ) %>%
  filter(.metric != "detection_prevalence")

supp.labs <-
  c(
    "Accuracy",
    "Kappa",
    "Sensitivity",
    "Specificity",
    "PPV",
    "NVP",
    "MCC",
    "J Index",
    "Balanced Accuracy",
    "Precision",
    "Recall",
    "F1 Score"
  )
names(supp.labs) <- unique(results_summarized$.metric)

results_summarized %>%
  ggplot(aes(x = s, y = .estimate, col = .metric)) +
  geom_point() +
  geom_line() +
  theme_light() +
  facet_wrap(. ~ .metric,
             scales = "free_y",
             labeller = labeller(.metric = supp.labs)) +
  labs(x = "Number of Seconds", y = "Estimate") +
  scale_x_continuous(breaks = seq(1, 21, 2)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "none")

rm(list = ls())

# CMA 
library(tidyverse)
library(purrr)
library(tidymodels)
tidymodels_prefer()
options(dplyr.summarise.inform=FALSE)

# functions for plotting/deriving CMA and unadjusted CIs 
# plots the cma and unadjusted grid cells for 3 lags: 15, 30, 45 centiseconds 
plot_w_points_cma <- 
  function(data_train, data_test, sub, raw_data){
    nzv_trans <- 
      recipe(ID ~ ., data = data_train) %>% 
      step_nzv(all_predictors())
    
    nzv_estimates <- prep(nzv_trans)
    
    nzv <- colnames(juice(nzv_estimates))
    dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv), -second)
    dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv), -second)
    
    train <- dat_nzv
    train$class <- ifelse(train$ID == sub, 1, 0)
    tmp <- train %>% dplyr::select(-c(ID)) 
    mod <- glm(class ~ ., data = tmp, family=binomial(link="logit"))
    
    
    summary <-
      mod %>% 
      tidy() %>%
      arrange(p.value) %>%
      filter(term != "(Intercept)") 
    
    terms <- summary$term
    # variance covariance matrix 
    # code to cma adjustment 
    v_hat <- vcov(mod)[terms, terms]
    a <- 1/(sqrt(diag(v_hat))) 
    A <- a %*% t(a)
    C <- v_hat*A
    
    q <- mvtnorm::qmvnorm(p = .95, corr = C)$quantile
    
    summary <- 
      summary %>% rowwise() %>% 
      mutate(
        lb_marg = estimate - (1.96*std.error),
        ub_marg = estimate + (1.96*std.error),
        lb_cma = estimate - (q*std.error),
        ub_cma = estimate +  (q*std.error),
        sig_marg = ifelse(between(0, lb_marg, ub_marg), 0, 1),
        sig_cma = ifelse(between(0, lb_cma, ub_cma), 0, 1),
      )
    summary <- 
      summary %>% rowwise() %>% 
      mutate(
        lag = str_sub(term, -3, -2),
        sig = sub('.', '', str_split(term, " ")[[1]][1]),
        lagsig = str_split(term, " ")[[1]][2]
      )
    
    if(nrow(summary %>% filter(sig_cma==1)) == 0){
      return(NULL)
    }
    else{
      train_times <-
        data_train %>%
        dplyr::select(ID, second) %>%
        filter(ID == sub)
      
      
      pts_dat <- 
        raw_data %>% inner_join(., train_times, by = c("ID" = "ID", "second" = "second")) %>%
        group_by(ID, second) %>% mutate(
          lag_1 = dplyr::lag(sig, 15),
          lag_2 = dplyr::lag(sig, 30),
          lag_3 = dplyr::lag(sig, 45)
        ) %>% dplyr::select(sig, lag_1, lag_2, lag_3, ID, second) %>% 
        pivot_longer(lag_1:lag_3) %>% filter(!is.na(value)) %>% 
        mutate(lag = case_when(
          name == "lag_1" ~ 15,
          name == "lag_2" ~ 30,
          name== "lag_3" ~ 45,
        ),
        xmin = value, 
        xmax = value,
        ymin = sig,
        ymax = sig) 
      
      labs <- c("u = 15", "u = 30", "u = 45")
      names(labs) <- c(15, 30, 45)
      g <- summary %>% 
        mutate(
          est = ifelse(sig_cma == 1, exp(estimate), NA)) %>% 
        mutate(
          xmin = as.numeric(sub(".*\\(", "", sub(",.*", "", lagsig))),
          xmax = as.numeric(str_sub(sub(".*,", "", lagsig), end = -2)),
          ymin =  as.numeric(sub(".*\\(", "", sub(",.*", "", sig))),
          ymax = as.numeric(str_sub(sub(".*,", "", sig), end = -2))
        ) %>% 
        ggplot(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+
        geom_rect(aes(fill = est)) +
        facet_wrap(.~lag, labeller = labeller(lag = labs))+
        scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                             midpoint = 1,na.value="white", name = "Odds Ratio")+
        labs(title = "",
             subtitle = "",
             x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), 
             y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))')
        )+
        scale_x_continuous(limits=c(0,3))+
        scale_y_continuous(limits=c(0,3))+
        theme_linedraw()+
        theme(strip.text = element_text(face = "italic", size = 15),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 15))
      
      g + geom_point(data = pts_dat, aes(x = value, y = sig), col = "black", 
                     alpha = 0.01, size = .4)
    }
  }

plot_w_points_unadj <- 
  function(data_train, data_test, sub, raw_data){
    nzv_trans <- 
      recipe(ID ~ ., data = data_train) %>% 
      step_nzv(all_predictors())
    
    nzv_estimates <- prep(nzv_trans)
    
    nzv <- colnames(juice(nzv_estimates))
    dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv), -second)
    dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv), -second)
    
    train <- dat_nzv
    train$class <- ifelse(train$ID == sub, 1, 0)
    tmp <- train %>% dplyr::select(-c(ID)) 
    mod <- glm(class ~ ., data = tmp, family=binomial(link="logit"))
    
    
    summary <-
      mod %>% 
      tidy() %>%
      arrange(p.value) %>%
      filter(term != "(Intercept)") 
    
    terms <- summary$term
    # variance covariance matrix 
    # code to cma adjustment 
    v_hat <- vcov(mod)[terms, terms]
    a <- 1/(sqrt(diag(v_hat))) 
    A <- a %*% t(a)
    C <- v_hat*A
    
    q <- mvtnorm::qmvnorm(p = .95, corr = C)$quantile
    
    summary <- 
      summary %>% rowwise() %>% 
      mutate(
        lb_marg = estimate - (1.96*std.error),
        ub_marg = estimate + (1.96*std.error),
        lb_cma = estimate - (q*std.error),
        ub_cma = estimate +  (q*std.error),
        sig_marg = ifelse(between(0, lb_marg, ub_marg), 0, 1),
        sig_cma = ifelse(between(0, lb_cma, ub_cma), 0, 1),
      )
    summary <- 
      summary %>% rowwise() %>% 
      mutate(
        lag = str_sub(term, -3, -2),
        sig = sub('.', '', str_split(term, " ")[[1]][1]),
        lagsig = str_split(term, " ")[[1]][2]
      )
    
    if(nrow(summary %>% filter(sig_marg==1)) == 0){
      return(NULL)
    }
    else{
      train_times <-
        data_train %>%
        dplyr::select(ID, second) %>%
        filter(ID == sub)
      
      
      pts_dat <- 
        raw_data %>% inner_join(., train_times, by = c("ID" = "ID", "second" = "second")) %>%
        group_by(ID, second) %>% mutate(
          lag_1 = dplyr::lag(sig, 15),
          lag_2 = dplyr::lag(sig, 30),
          lag_3 = dplyr::lag(sig, 45)
        ) %>% dplyr::select(sig, lag_1, lag_2, lag_3, ID, second) %>% 
        pivot_longer(lag_1:lag_3) %>% filter(!is.na(value)) %>% 
        mutate(lag = case_when(
          name == "lag_1" ~ 15,
          name == "lag_2" ~ 30,
          name== "lag_3" ~ 45,
        ),
        xmin = value, 
        xmax = value,
        ymin = sig,
        ymax = sig) 
      
      labs <- c("u = 15", "u = 30", "u = 45")
      names(labs) <- c(15, 30, 45)
      g <- summary %>% 
        mutate(
          est = ifelse(sig_marg == 1, exp(estimate), NA)) %>% 
        mutate(
          xmin = as.numeric(sub(".*\\(", "", sub(",.*", "", lagsig))),
          xmax = as.numeric(str_sub(sub(".*,", "", lagsig), end = -2)),
          ymin =  as.numeric(sub(".*\\(", "", sub(",.*", "", sig))),
          ymax = as.numeric(str_sub(sub(".*,", "", sig), end = -2))
        ) %>% 
        ggplot(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+
        geom_rect(aes(fill = est)) +
        facet_wrap(.~lag, labeller = labeller(lag = labs))+
        scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                             midpoint = 1,na.value="white", name = "Odds Ratio")+
        labs(title = "",
             subtitle = "",
             x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), 
             y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))')
        )+
        scale_x_continuous(limits=c(0,3))+
        scale_y_continuous(limits=c(0,3))+
        theme_linedraw()+
        theme(strip.text = element_text(face = "italic", size = 15),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 15))
      
      g + geom_point(data = pts_dat, aes(x = value, y = sig), col = "black", 
                     alpha = 0.01, size = .4)
    }
  }

## for IU data 
grid_data_lw_IU <- readRDS(here::here("data/grid_data_lw_IU.rds"))
df_all_IU <- read_csv(here::here("data/df_all_IU.csv"), 
                      col_types = cols(...1 = col_skip())) %>%
  rename(ID = ID2,
         signal_lwrist = signal_lw,
         signal_lhip = signal_lh,
         signal_lankle = signal_la,
         signal_rankle = signal_ra)

# function to get training and testing data 
get_train_test <- function(pct, data){
  data_split <- split(data, f = data$ID)
  
  samp <- function(pct, n, ind) {
    set.seed(ind)
    sample(n, floor(pct * n), replace = F)
  }
  ids <- unique(data$ID)
  # number of rows for each individual 
  rows <- lapply(data_split, nrow) %>% unlist()
  # get random 75% of seconds for training 
  train_indices <- map2(pct = pct,
                        .x = rows,
                        .y = ids,
                        .f = samp)
  
  getrows <- function(data, rows) {
    data[rows, ]
  }
  getrows_test <- function(data, rows) {
    data[-rows, ]
  }
  
  data_train <-
    map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
  data_test <-
    map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)
  return(list(train = data_train, test = data_test))
}

data_train_IU <- get_train_test(pct = .75, data = grid_data_lw_IU)$train
data_test_IU <- get_train_test(pct = .75, data = grid_data_lw_IU)$test


plot_w_points_cma(data_train_IU, data_test_IU, sub = 3, raw_data = df_all_IU %>%
                    rename(sig = signal_lwrist))
plot_w_points_unadj(data_train_IU, data_test_IU, sub = 3, raw_data = df_all_IU %>%
                      rename(sig = signal_lwrist))


grid_data_rw_zju_s1 <- readRDS(here::here("data/grid_data_rw_zju_s1.rds")) %>%
  mutate(ID = ID - 22)
df_all_zju <- read_csv(here::here("data/df_all_zju.csv"), 
                       col_types = cols(...1 = col_skip())) %>%
  filter(session == "session_1") %>%
  mutate(ID = ID - 22) 

data_train_zju <- get_train_test(pct = .75, data = grid_data_rw_zju_s1)$train
data_test_zju <- get_train_test(pct = .75, data = grid_data_rw_zju_s1)$test

plot_w_points_cma(data_train_zju, data_test_zju, sub = 143, raw_data = df_all_zju %>%
                    rename(sig = signal_rwrist))
plot_w_points_unadj(data_train_zju, data_test_zju, sub = 143, raw_data = df_all_zju %>%
                      rename(sig = signal_rwrist))

# function to compare the number of sig grid cells - marginal and CMA 
compare_cma_unadj <-
  function(data_train, data_test, sub){
    nzv_trans <- 
      recipe(ID ~ ., data = data_train) %>% 
      step_nzv(all_predictors())
    
    nzv_estimates <- prep(nzv_trans)
    
    nzv <- colnames(juice(nzv_estimates))
    dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv), -second)
    dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv), -second)
    
    train <- dat_nzv
    train$class <- ifelse(train$ID == sub, 1, 0)
    tmp <- train %>% dplyr::select(-c(ID)) 
    mod <- glm(class ~ ., data = tmp, family=binomial(link="logit"))
    
    
    summary <-
      mod %>% 
      tidy() %>%
      arrange(p.value) %>%
      filter(term != "(Intercept)") 
    
    terms <- summary$term
    # variance covariance matrix 
    # code to cma adjustment 
    v_hat <- vcov(mod)[terms, terms]
    a <- 1/(sqrt(diag(v_hat))) 
    A <- a %*% t(a)
    C <- v_hat*A
    
    q <- mvtnorm::qmvnorm(p = .95, corr = C)$quantile
    
    summary <- 
      summary %>% rowwise() %>% 
      mutate(
        lb_marg = estimate - (1.96*std.error),
        ub_marg = estimate + (1.96*std.error),
        lb_cma = estimate - (q*std.error),
        ub_cma = estimate +  (q*std.error),
        sig_marg = ifelse(between(0, lb_marg, ub_marg), 0, 1),
        sig_cma = ifelse(between(0, lb_cma, ub_cma), 0, 1),
      )
    summary <- 
      summary %>% rowwise() %>% 
      mutate(
        lag = str_sub(term, -3, -2),
        sig = sub('.', '', str_split(term, " ")[[1]][1]),
        lagsig = str_split(term, " ")[[1]][2]
      )
    
    return(data.frame(sub = sub,
                      n_marg = sum(summary$sig_marg),
                      n_cma = sum(summary$sig_cma),
                      n_total = nrow(summary)))
  }

comparison <-
  map(.x = seq(1, 32, 1),
      .f = compare_cma_unadj,
      data_train = data_train_IU,
      data_test = data_test_IU) %>%
  list_rbind()


comparison_zju <-
  map(.x = seq(1, 153, 1),
      .f = compare_cma_unadj,
      data_train = data_train_zju,
      data_test = data_test_zju) %>%
  list_rbind()

rm(list = ls())
# IU machine learning
library(tidymodels)
library(tune)
library(baguette)
library(discrim)
library(finetune)
library(readr)
tidymodels_prefer()

data <- readRDS(here::here("data/grid_data_lw_IU.rds")) %>%
  janitor::clean_names() %>%
  rename("ID" = "id")
data_split <- split(data, f = factor(data$ID))

# function to get random sample
samp <- function(pct, n, ind) {
  set.seed(ind)
  sample(n, floor(pct * n), replace = F)
}
ids <- unique(data$ID)
rows <- lapply(data_split, nrow) %>% unlist()

train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

getrows <- function(data, rows) {
  data[rows,]
}
getrows_test <- function(data, rows) {
  data[-rows,]
}
data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


get_original_indices <- function(subject, list) {
  cbind(rep(subject, length(list[[subject]])), list[[subject]]) %>%
    data.frame() %>%
    rename(sub = X1,
           sec = X2)
}

orig_indices_key <-
  map_dfr(.x = ids,
          .f = get_original_indices,
          list = train_indices)

train_inds <-
  data %>%
  full_join(orig_indices_key %>% mutate(train = 1),
            by = c("ID" = "sub",
                   "second" = "sec")) %>%
  ungroup() %>%
  mutate(row = row_number())

true_train <- train_inds$row[!is.na(train_inds$train)]
# create folds for cross validation


# initial split data for eventually testing
data_split <- initial_split(data, prop = 0.75)
data_split$in_id <- true_train

# check
#all_equal(data_train, training(data_split2))
#all_equal(data_test, testing(data_split2))

data_split$data <- data %>% dplyr::select(-second)
data_train <- data_train %>% dplyr::select(-second)
data_test <- data_test %>% dplyr::select(-second)
data_folds <- vfold_cv(data_train, v = 5)
## model specifications

# flexible discriminant model from multivariate adaptive regression splines
discrim_flex_spec <-
  discrim_flexible(num_terms = tune(), prod_degree = tune()) %>%
  set_engine("earth") %>%
  set_mode("classification")

# neural net
nnet_spec <-
  mlp(hidden_units = tune(),
      penalty = tune(),
      epochs = tune()) %>%
  set_engine("nnet", MaxNWts = 2600) %>%
  set_mode("classification")

# radial basis function support vector machine
svm_r_spec <-
  svm_rbf(cost = tune(),
          rbf_sigma = tune(),
          margin = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# polynomial basis function support vector machine
svm_p_spec <-
  svm_poly(
    cost = tune(),
    degree = tune(),
    scale_factor = tune(),
    margin = tune()
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# bayesian additive regression tree
bart_spec <-
  parsnip::bart(
    trees = tune(),
    prior_terminal_node_coef = tune(),
    prior_terminal_node_expo = tune()
  ) %>%
  set_engine("dbarts") %>%
  set_mode("classification")

# boosted tree
boost_spec <-
  boost_tree(
    min_n = tune(),
    mtry = tune(),
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    sample_size = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# multivariate adaptive regression spline
mars_spec <-
  mars(num_terms = tune(), prod_degree = tune()) %>%
  set_engine("earth") %>%
  set_mode("classification")

# k nearest neighbors
knn_spec <-
  nearest_neighbor(neighbors = tune()) %>%
  set_engine("kknn") %>%
  set_mode("classification")

# logistic regression
lreg_spec <-
  logistic_reg(penalty = 0, mixture = 0) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# penalized logistic regression
plreg_spec <-
  logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# naive bayes
bayes_spec <-
  naive_Bayes(smoothness = tune()) %>%
  set_engine("naivebayes") %>%
  set_mode("classification")

# random forest
rf_spec <-
  rand_forest(mtry = tune(),
              trees = 1000,
              min_n = tune()) %>%
  set_engine("ranger") %>%
  set_mode("classification")

# function that finds and fits best model for each subject
get_predictions_models <- function(subject) {
  tmp_train <-
    data_train %>% mutate(ID = factor(ifelse(ID == subject, 1, 0)))
  tmp_folds <- vfold_cv(tmp_train, v = 5)
  
  tmp_data <- data %>% mutate(ID = factor(ifelse(ID == subject, 1, 0)))
  
  data_split_tmp <- data_split
  data_split_tmp$data <- tmp_data
  
  
  # recipes
  normalized_transform_rec <-
    recipe(ID ~ ., data = tmp_train) %>%
    step_nzv(all_predictors()) %>%
    step_YeoJohnson(all_predictors()) %>%
    step_normalize(all_predictors())
  
  no_rec <-
    recipe(ID ~ ., data = tmp_train) %>%
    step_nzv(all_predictors())
  
  # workflows
  proc  <- workflow_set(
    preproc = list(all = normalized_transform_rec),
    models = list(
      mlp = nnet_spec,
      svm_rf = svm_r_spec,
      svm_p = svm_p_spec,
      knn = knn_spec
    )
  )
  
  none <- workflow_set(
    preproc = list(none = no_rec),
    models = list(
      mars  = mars_spec,
      bart = bart_spec,
      boost = boost_spec,
      flexdiscrim = discrim_flex_spec,
      rf = rf_spec,
      naivebayes = bayes_spec,
      lreg = lreg_spec,
      plreg = plreg_spec
    )
  )
  
  all_workflows <-
    bind_rows(proc, none)
  
  
  race_ctrl <-
    control_race(
      save_pred = TRUE,
      parallel_over = "everything",
      save_workflow = TRUE
    )
  
  results_all <-
    all_workflows %>%
    workflow_map(
      "tune_race_anova",
      seed = 1503,
      resamples = tmp_folds,
      grid = 25,
      control = race_ctrl
    )
  
  results <-
    results_all %>%
    filter(lengths(result) > 1) %>%
    workflowsets::rank_results(select_best = TRUE, rank_metric = "roc_auc") %>%
    select(wflow_id, .metric, race = mean, config_race = .config) %>%
    filter(.metric == "roc_auc")
  
  best_id <-
    results_all %>%
    filter(lengths(result) > 1) %>%
    workflowsets::rank_results(select_best = TRUE, rank_metric = "roc_auc") %>%
    select(wflow_id, .metric, race = mean, config_race = .config) %>%
    filter(.metric == "roc_auc") %>%
    mutate(rank = rank(-1 * race)) %>% filter(rank == 1) %>%
    dplyr::select(wflow_id) %>% unlist() %>% unname()
  
  best_results <-
    results_all %>%
    extract_workflow_set_result(id = best_id) %>%
    select_best(metric = "roc_auc")
  
  test_results <-
    results_all %>%
    extract_workflow(id = best_id) %>%
    finalize_workflow(best_results) %>%
    last_fit(split = data_split_tmp)
  
  preds <- test_results %>%
    collect_predictions() %>%
    mutate(model = best_id,
           subject = subject)
  return(list(preds = preds, results = results_all))
}

subs <-
  data %>%
  dplyr::select(ID) %>%
  unlist() %>%
  unname() %>%
  unique() %>%
  as.character()

# to run entire thing at once
# can break this up because it takes a while to run

all_predictions_models <-
  map(.x = subs, .f = get_predictions_models)


ranked_results <- function(id, list) {
  list[[id]][2]$results %>%
    filter(lengths(result) > 1) %>%
    rank_results(rank_metric = "roc_auc", select_best = TRUE) %>%
    mutate(subject = id)
}

all_results <-
  map_dfr(.x = subs,
          .f = ranked_results,
          list = all_predictions_models)


get_preds <- function(id, list) {
  list[[id]][1]$preds %>%
    arrange(.row) %>%
    dplyr::select(.pred_1)
}

all_preds <-
  map_dfc(.x = subs,
          .f = get_preds,
          list = all_predictions_models)

if (!dir.exists("predictions"))
  dir.create(here::here("predictions"))


colnames(all_preds) <- paste0("x", seq(1, 32, 1))

row_sums <- rowSums(all_preds)

# normalize and add "true subject column"
all_preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x32, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = data_test$ID)


saveRDS(all_preds,
        here::here("predictions/IU_ml_predictions.rds"))

get_summarized_predictions <- function(predictions, long = FALSE) {
  if (long == T) {
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      mutate(correct = ifelse(true_subject == model, 1, 0))
  }
  else{
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      group_by(true_subject) %>%
      summarize(
        maxprob = first(max(mean_pred)),
        predicted_sub = first(model[mean_pred == maxprob]),
        probsubj = first(mean_pred[true_subject == model])
      ) %>%
      mutate(correct = ifelse(as.numeric(predicted_sub) == true_subject, 1, 0))
  }
}

get_summarized_predictions(all_preds) %>% print(n = Inf)


saveRDS(all_results,
        here::here("predictions/IU_ml_model_results.rds"))
rm(list = ls())

# zju machine learning 
library(tidymodels)
library(tune)
library(baguette)
library(discrim)
library(finetune)
library(readr)
library(doMC)
doMC::registerDoMC(cores = max(1, parallelly::availableCores() - 1))
tidymodels_prefer()


# session 1 training and testing 
# change column names to avoid bugs with some of the models 
data <- readRDS(here::here("data/grid_data_rw_zju_s1.rds")) %>% 
  mutate(
    ID = ID - 22
  ) %>%
  janitor::clean_names() %>%
  rename(
    ID = id) 


data_split <- split(data, f = factor(data$ID))

# function to get random sample 
samp <- function(pct, n, ind) {
  set.seed(ind)
  sample(n, floor(pct * n), replace = F)
}
ids <- unique(data$ID)
rows <- lapply(data_split, nrow) %>% unlist()

train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

getrows <- function(data, rows) {
  data[rows, ]
}
getrows_test <- function(data, rows) {
  data[-rows, ]
}
data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


get_original_indices <- function(subject, list){
  cbind(rep(subject, length(list[[subject]])), list[[subject]]) %>% 
    data.frame() %>%
    rename(
      sub = X1, 
      sec = X2
    )
}

orig_indices_key <- 
  map_dfr(.x = ids,
          .f = get_original_indices,
          list = train_indices)

train_inds <-
  data %>%
  full_join(orig_indices_key %>% mutate(train = 1), by = c("ID" = "sub",
                                                           "second" = "sec")) %>%
  ungroup() %>% 
  mutate(
    row = row_number()
  )

true_train <- train_inds$row[!is.na(train_inds$train)]
# create folds for cross validation 


# initial split data for eventually testing 
data_split <- initial_split(data, prop = 0.75)
data_split$in_id <- true_train

# check 
#all_equal(data_train, training(data_split2))
#all_equal(data_test, testing(data_split2))

data_split$data <- data %>% dplyr::select(-second)
data_train <- data_train %>% dplyr::select(-second)
data_test <- data_test %>% dplyr::select(-second)
data_folds <- vfold_cv(data_train, v = 5)
## model specifications 

# flexible discriminant model from multivariate adaptive regression splines 
discrim_flex_spec <-
  discrim_flexible(num_terms = tune(), prod_degree = tune()) %>%
  set_engine("earth") %>%
  set_mode("classification")

# neural net 
nnet_spec <- 
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>% 
  set_engine("nnet", MaxNWts = 2600) %>% 
  set_mode("classification")

# radial basis function support vector machine 
svm_r_spec <- 
  svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")

# polynomial basis function support vector machine 
svm_p_spec <- 
  svm_poly(cost = tune(), degree = tune(), scale_factor = tune(), margin = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")

# bayesian additive regression tree 
bart_spec <- 
  parsnip::bart(trees=tune(), prior_terminal_node_coef = tune(), prior_terminal_node_expo = tune()) %>% 
  set_engine("dbarts") %>%
  set_mode("classification")

# boosted tree 
boost_spec <- 
  boost_tree(min_n = tune(), mtry = tune(), trees = tune(), tree_depth = tune(),
             learn_rate = tune(), sample_size = tune(), loss_reduction = tune()) %>% 
  set_engine("xgboost") %>%
  set_mode("classification")

# multivariate adaptive regression spline 
mars_spec <- 
  mars(num_terms = tune(), prod_degree = tune()) %>%
  set_engine("earth") %>%
  set_mode("classification")

# k nearest neighbors 
knn_spec <- 
  nearest_neighbor(neighbors = tune()) %>%
  set_engine("kknn") %>% 
  set_mode("classification")

# logistic regression 
lreg_spec <-
  logistic_reg(penalty = 0, mixture = 0) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# penalized logistic regression 
plreg_spec <-
  logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# naive bayes 
bayes_spec <-
  naive_Bayes(smoothness = tune()) %>%
  set_engine("naivebayes") %>%
  set_mode("classification")

# random forest 
rf_spec <- 
  rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")

# function that finds and fits best model for each subject 
get_predictions_models <- function(subject){
  
  race_ctrl <-
    control_race(
      save_pred = F,
      parallel_over = "everything",
      save_workflow = F
    )
  
  
  
  tmp_train <- 
    data_train %>% mutate(
      ID = factor(ifelse(ID == subject, 1, 0))
    )
  tmp_folds <- vfold_cv(tmp_train, v = 5)
  
  tmp_data <- data %>% mutate(
    ID = factor(ifelse(ID == subject, 1, 0))
  )
  
  data_split_tmp <- data_split
  data_split_tmp$data <- tmp_data 
  
  
  # recipes 
  normalized_transform_rec <- 
    recipe(ID ~ ., data = tmp_train) %>% 
    step_nzv(all_predictors()) %>% 
    step_YeoJohnson(all_predictors()) %>%
    step_normalize(all_predictors())
  
  no_rec <- 
    recipe(ID ~ ., data = tmp_train) %>% 
    step_nzv(all_predictors()) 
  
  # workflows 
  proc  <- workflow_set(
    preproc = list(all = normalized_transform_rec),      
    models = list(mlp = nnet_spec, svm_rf = svm_r_spec, svm_p = svm_p_spec, knn = knn_spec)
  )
  
  none <- workflow_set(
    preproc = list(none = no_rec),
    models = list(mars  = mars_spec, bart = bart_spec, boost = boost_spec, flexdiscrim = discrim_flex_spec,
                  rf = rf_spec, naivebayes = bayes_spec, lreg = lreg_spec, plreg = plreg_spec)
  )
  
  all_workflows <- 
    bind_rows(proc, none) 
  
  results_all <-
    all_workflows %>% 
    workflow_map(
      "tune_race_anova",
      seed = 1503,
      resamples = tmp_folds,
      grid = 25,
      control = race_ctrl
    )
  
  results <- 
    results_all %>% 
    filter(lengths(result) > 1) %>% 
    workflowsets::rank_results(select_best = TRUE, rank_metric = "roc_auc") %>% 
    select(wflow_id, .metric, race = mean, config_race = .config) %>%
    filter(.metric == "roc_auc")
  
  best_id <- 
    results_all %>% 
    filter(lengths(result) > 1) %>% 
    workflowsets::rank_results(select_best = TRUE, rank_metric = "roc_auc") %>% 
    select(wflow_id, .metric, race = mean, config_race = .config) %>%
    filter(.metric == "roc_auc") %>%
    mutate(
      rank = rank(-1*race)
    ) %>% filter(rank == 1) %>%
    dplyr::select(wflow_id) %>% unlist() %>% unname()
  
  best_results <- 
    results_all %>%
    extract_workflow_set_result(id = best_id) %>%
    select_best(metric = "roc_auc")
  
  test_results <- 
    results_all %>%
    extract_workflow(id = best_id) %>%
    finalize_workflow(best_results) %>%
    last_fit(split = data_split_tmp)
  
  preds <- test_results %>% 
    collect_predictions() %>% 
    mutate(
      model = best_id,
      subject = subject
    )
  return(list(preds = preds, results = results_all))
}

subs <- 
  data %>% 
  dplyr::select(ID) %>%
  unlist() %>% 
  unname() %>% 
  unique() %>%
  as.character() 

all_predictions_models <- 
  map(.x = subs, 
      .f = get_predictions_models,
      .progress = T)

ranked_results <- function(id, list) {
  list[[id]][2]$results %>%
    filter(lengths(result) > 1) %>%
    rank_results(rank_metric = "roc_auc", select_best = TRUE) %>%
    mutate(subject = id)
}

all_results <-
  map_dfr(.x = subs,
          .f = ranked_results,
          list = all_predictions_models)


get_preds <- function(id, list) {
  list[[id]][1]$preds %>%
    arrange(.row) %>%
    dplyr::select(.pred_1)
}

all_preds <-
  map_dfc(.x = subs,
          .f = get_preds,
          list = all_predictions_models)

if (!dir.exists("predictions"))
  dir.create(here::here("predictions"))


colnames(all_preds) <- paste0("x", seq(1, 153, 1))

row_sums <- rowSums(all_preds)

# normalize and add "true subject column"
all_preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = data_test$ID)


saveRDS(all_preds,
        here::here("predictions/zjus1_ml_predictions.rds"))

get_summarized_predictions <- function(predictions, long = FALSE) {
  if (long == T) {
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      mutate(correct = ifelse(true_subject == model, 1, 0))
  }
  else{
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      group_by(true_subject) %>%
      summarize(
        maxprob = first(max(mean_pred)),
        predicted_sub = first(model[mean_pred == maxprob]),
        probsubj = first(mean_pred[true_subject == model])
      ) %>%
      mutate(correct = ifelse(as.numeric(predicted_sub) == true_subject, 1, 0))
  }
}

get_summarized_predictions(all_preds) %>% print(n = Inf)

saveRDS(all_results,
        here::here("predictions/zjus1_ml_model_results.rds"))

# session 1 training, session 2 testing 
data_s1 <- readRDS(here::here("data/grid_data_rw_zju_s1.rds")) %>% 
  mutate(
    ID = ID - 22
  ) %>%
  janitor::clean_names() %>%
  rename(
    ID = id
  )

data_s2 <- readRDS(here::here("data/grid_data_rw_zju_s2.rds")) %>% 
  mutate(
    ID = ID - 22
  ) %>%
  janitor::clean_names() %>%
  rename(
    ID = id
  )


## first 10024 rows are s1
data <- bind_rows(data_s1, data_s2)

# initial split data for eventually testing 
data_split <- initial_split(data, prop = 0.75)
data_split$in_id <- 1:10024

# check 
#all_equal(data_train, training(data_split2))
#all_equal(data_test, testing(data_split2))

data_split$data <- data %>% dplyr::select(-second)
data_train <- data_s1 %>% dplyr::select(-second)
data_test <- data_s2 %>% dplyr::select(-second)
data_folds <- vfold_cv(data_train, v = 5)
## model specifications 

# flexible discriminant model from multivariate adaptive regression splines 
discrim_flex_spec <-
  discrim_flexible(num_terms = tune(), prod_degree = tune()) %>%
  set_engine("earth") %>%
  set_mode("classification")

# neural net 
nnet_spec <- 
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>% 
  set_engine("nnet", MaxNWts = 2600) %>% 
  set_mode("classification")

# radial basis function support vector machine 
svm_r_spec <- 
  svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")

# polynomial basis function support vector machine 
svm_p_spec <- 
  svm_poly(cost = tune(), degree = tune(), scale_factor = tune(), margin = tune()) %>% 
  set_engine("kernlab") %>% 
  set_mode("classification")

# bayesian additive regression tree 
bart_spec <- 
  parsnip::bart(trees=tune(), prior_terminal_node_coef = tune(), prior_terminal_node_expo = tune()) %>% 
  set_engine("dbarts") %>%
  set_mode("classification")

# boosted tree 
boost_spec <- 
  boost_tree(min_n = tune(), mtry = tune(), trees = tune(), tree_depth = tune(),
             learn_rate = tune(), sample_size = tune(), loss_reduction = tune()) %>% 
  set_engine("xgboost") %>%
  set_mode("classification")

# multivariate adaptive regression spline 
mars_spec <- 
  mars(num_terms = tune(), prod_degree = tune()) %>%
  set_engine("earth") %>%
  set_mode("classification")

# k nearest neighbors 
knn_spec <- 
  nearest_neighbor(neighbors = tune()) %>%
  set_engine("kknn") %>% 
  set_mode("classification")

# logistic regression 
lreg_spec <-
  logistic_reg(penalty = 0, mixture = 0) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# penalized logistic regression 
plreg_spec <-
  logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# naive bayes 
bayes_spec <-
  naive_Bayes(smoothness = tune()) %>%
  set_engine("naivebayes") %>%
  set_mode("classification")

# random forest 
rf_spec <- 
  rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")

# function that finds and fits best model for each subject 
get_predictions_models <- function(subject){
  
  race_ctrl <-
    control_race(
      save_pred = F,
      parallel_over = "everything",
      save_workflow = TRUE
    )
  tmp_train <- 
    data_train %>% mutate(
      ID = factor(ifelse(ID == subject, 1, 0))
    )
  tmp_folds <- vfold_cv(tmp_train, v = 5)
  
  tmp_data <- data %>% mutate(
    ID = factor(ifelse(ID == subject, 1, 0))
  )
  
  data_split_tmp <- data_split
  data_split_tmp$data <- tmp_data 
  
  
  # recipes 
  normalized_transform_rec <- 
    recipe(ID ~ ., data = tmp_train) %>% 
    step_nzv(all_predictors()) %>% 
    step_YeoJohnson(all_predictors()) %>%
    step_normalize(all_predictors())
  
  no_rec <- 
    recipe(ID ~ ., data = tmp_train) %>% 
    step_nzv(all_predictors()) 
  
  # workflows 
  proc  <- workflow_set(
    preproc = list(all = normalized_transform_rec),      
    models = list(mlp = nnet_spec, svm_rf = svm_r_spec, svm_p = svm_p_spec, knn = knn_spec)
  )
  
  none <- workflow_set(
    preproc = list(none = no_rec),
    models = list(mars  = mars_spec, bart = bart_spec, boost = boost_spec, flexdiscrim = discrim_flex_spec,
                  rf = rf_spec, naivebayes = bayes_spec, lreg = lreg_spec, plreg = plreg_spec)
  )
  
  all_workflows <- 
    bind_rows(proc, none) 
  
  results_all <-
    all_workflows %>% 
    workflow_map(
      "tune_race_anova",
      seed = 1503,
      resamples = tmp_folds,
      grid = 25,
      control = race_ctrl
    )
  
  results <- 
    results_all %>% 
    filter(lengths(result) > 1) %>% 
    workflowsets::rank_results(select_best = TRUE, rank_metric = "roc_auc") %>% 
    select(wflow_id, .metric, race = mean, config_race = .config) %>%
    filter(.metric == "roc_auc")
  
  best_id <- 
    results_all %>% 
    filter(lengths(result) > 1) %>% 
    workflowsets::rank_results(select_best = TRUE, rank_metric = "roc_auc") %>% 
    select(wflow_id, .metric, race = mean, config_race = .config) %>%
    filter(.metric == "roc_auc") %>%
    mutate(
      rank = rank(-1*race)
    ) %>% filter(rank == 1) %>%
    dplyr::select(wflow_id) %>% unlist() %>% unname()
  
  best_results <- 
    results_all %>%
    extract_workflow_set_result(id = best_id) %>%
    select_best(metric = "roc_auc")
  
  test_results <- 
    results_all %>%
    extract_workflow(id = best_id) %>%
    finalize_workflow(best_results) %>%
    last_fit(split = data_split_tmp)
  
  preds <- test_results %>% 
    collect_predictions() %>% 
    mutate(
      model = best_id,
      subject = subject
    )
  return(list(preds = preds, results = results_all))
}

subs <- 
  data_s1 %>% 
  dplyr::select(ID) %>%
  unlist() %>% 
  unname() %>% 
  unique() %>%
  as.character() 

all_predictions_models <- 
  map(.x = subs, 
      .f = get_predictions_models,
      .progress = T)

ranked_results <- function(id, list) {
  list[[id]][2]$results %>%
    filter(lengths(result) > 1) %>%
    rank_results(rank_metric = "roc_auc", select_best = TRUE) %>%
    mutate(subject = id)
}

all_results <-
  map_dfr(.x = subs,
          .f = ranked_results,
          list = all_predictions_models)


get_preds <- function(id, list) {
  list[[id]][1]$preds %>%
    arrange(.row) %>%
    dplyr::select(.pred_1)
}

all_preds <-
  map_dfc(.x = subs,
          .f = get_preds,
          list = all_predictions_models)


colnames(all_preds) <- paste0("x", seq(1, 153, 1))

row_sums <- rowSums(all_preds)

# normalize and add "true subject column"
all_preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = data_test$ID)


saveRDS(all_preds,
        here::here("predictions/zjus1s2_ml_predictions.rds"))


get_summarized_predictions(all_preds) %>% print(n = Inf)

saveRDS(all_results,
        here::here("predictions/zjus1s2_ml_model_results.rds"))

rm(list = ls())

# iu functional
library(tidyverse)
library(mgcv)
library(readr)
library(magrittr)
expit <- function(x)
  1 / (1 + exp(-x))
`%notin%` <- Negate(`%in%`)


df_all_IU <- read_csv(here::here("data/df_all_IU.csv"),
                      col_types = cols(...1 = col_skip())) %>%
  rename(
    ID = ID2,
    signal_lwrist = signal_lw,
    signal_lhip = signal_lh,
    signal_lankle = signal_la,
    signal_rankle = signal_ra
  )

# tlen is the length of intervals here, we choose 100 centiseconds
get_functional_data <- function(df, tlen) {
  uid <- unique(df$ID)
  nid <- length(uid)
  
  nk <- (tlen * 100) * (tlen * 100 - 1) / 2
  ji <- df %>%
    group_by(ID) %>%
    summarize(maxtime = max(second)) %>%
    dplyr::select(maxtime) %>%
    unlist() %>%
    unname()
  
  N <- sum(ji)
  smat_w <-
    umat_w <- dmat <- matrix(NA, ncol = nk, nrow = N) # matrix of NAs
  id_vec <- time_g_vec <- rep(NA, N)
  
  inx_samp <- lapply(1:(100 * tlen - 1), function(x) {
    ret <- data.frame("u" = x, "s" = (x + 1):(100 * tlen))
    ret
  })
  inx_samp <- bind_rows(inx_samp)
  inx <- 1
  
  for (i in seq_along(uid)) {
    df_i <- filter(df, ID == i)
    j_i    <- unique(df_i$second)
    for (j in seq_along(j_i)) {
      Y_ij <- df_i$signal_lwrist[df_i$second == j_i[j]]
      umat_w[inx, ] <- Y_ij[inx_samp$u]
      smat_w[inx, ] <- Y_ij[inx_samp$s]
      dmat[inx, ] <- inx_samp$s - inx_samp$u
      time_g_vec[inx] <- ji[j]
      id_vec[inx] <- df_i$ID[1]
      inx <- inx + 1
    }
  }
  
  df_fit <- data.frame(
    "ID" = id_vec,
    "umat" = I(umat_w),
    "smat" = I(smat_w),
    "dmat" = I(dmat),
    "lmat" = I(matrix(1 / nk, ncol = nk, nrow = N))
  )
  rm(list = c("smat_w", "umat_w", "dmat"))
  gc()
  
  
  df_fit <- df_fit %>%
    group_by(ID) %>%
    dplyr::mutate(J = 1:n()) %>%
    ungroup()
  
  df_fit
}

df_fit_IU <- get_functional_data(df_all_IU, tlen = 1)
data_split <- split(df_fit_IU, f = df_fit_IU$ID)

# function to sample percentage
samp <- function(pct, n, ind) {
  set.seed(ind)
  sample(n, floor(pct * n), replace = F)
}
ids <- unique(df_fit_IU$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

getrows <- function(data, rows) {
  subset(data, J %in% rows)
}
getrows_test <- function(data, rows) {
  subset(data, J %notin% rows)
}

data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


fit_functional_model <- function(train, test) {
  uid        <- unique(train$ID)
  nid        <- length(uid)
  
  lpmat <-
    matrix(NA, ncol = nid, nrow = nrow(test)) # to store predictions
  for (i in 1:nid) {
    train$Y <- as.numeric(train$ID == uid[i])
    fit_i <-
      gam(
        Y ~ te(umat, smat, dmat, by = lmat),
        method = "REML",
        data = train,
        family = quasibinomial()
      )
    lp_i <- predict.gam(fit_i, newdata = test, type = "link")
    lpmat[, i] <- unname(lp_i)
    rm(list = c("fit_i", "lp_i"))
    train$Y <- NULL
    print(i)
  }
  lpmat
}


## NOTE: this will likely not run on local machine, need to run on computing cluster
preds <- fit_functional_model(data_train, data_test)

preds %<>%
  data.frame() %>%
  janitor::clean_names() %>%
  expit()

row_sums <- rowSums(preds)

# normalize and add "true subject column"

preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x32, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = data_test$ID)

if (!dir.exists("predictions")) {
  dir.create(here::here("predictions"))
}

saveRDS(preds, here::here("predictions/IU_func_predictions.rds"))

# function to get prediction summaries
get_summarized_predictions <- function(predictions, long = FALSE) {
  if (long == T) {
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      mutate(correct = ifelse(true_subject == model, 1, 0))
  }
  else{
    predictions %>%
      group_by(true_subject) %>%
      mutate(sec = row_number()) %>%
      pivot_longer(cols = -c("true_subject", "sec")) %>%
      mutate(model = as.numeric(sub(".*x", "", name))) %>%
      rename(pred = value) %>%
      ungroup() %>%
      group_by(true_subject, model) %>%
      summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
      group_by(true_subject) %>%
      summarize(
        maxprob = first(max(mean_pred)),
        predicted_sub = first(model[mean_pred == maxprob]),
        probsubj = first(mean_pred[true_subject == model])
      ) %>%
      mutate(correct = ifelse(as.numeric(predicted_sub) == true_subject, 1, 0))
  }
}

# accuracy
get_summarized_predictions(preds) %>% print(n = Inf)

# perfectly predicted

get_accuracies <- function(predictions, seconds) {
  tmp <-
    predictions %>%
    group_by(true_subject) %>%
    mutate(sec = floor(row_number() / seconds)) %>%
    pivot_longer(cols = -c("true_subject", "sec")) %>%
    mutate(model = as.numeric(sub(".*x", "", name))) %>%
    rename(pred = value) %>%
    ungroup() %>%
    group_by(true_subject, model, sec) %>%
    summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
    group_by(true_subject, sec) %>%
    mutate(rank = rank(-mean_pred)) %>%
    filter(model == true_subject) %>%
    mutate(rank1 = ifelse(rank == 1, 1, 0),
           rank5 = ifelse(rank <= 5, 1, 0))
  n = nrow(tmp)
  r1 = sum(tmp$rank1) / n
  r5 = sum(tmp$rank5) / n
  return(data.frame(
    metric = c("rank1", "rank5"),
    value = c(r1, r5),
    s = c(seconds, seconds)
  ))
  
}


secs <- c(1, 2, 5, 10, 15, 20, 25, 30, 35, 50, 60, 70, 80, 90, 100)
stats_fnl <-
  map_dfr(.x = secs,
          .f = get_accuracies,
          predictions = preds) %>%
  mutate(model = "functional")


rm(list = ls())

# zju functional 
library(tidyverse)
library(mgcv)
library(readr)
library(purrr)
`%notin%` <- Negate(`%in%`)
expit <- function(x)
  1 / (1 + exp(-x))


# train and test on session 1
df_all_zju_s1 <- read_csv(here::here("data/df_all_zju.csv"),
                          col_types = cols(...1 = col_skip())) %>%
  filter(session == "session_1") %>%
  mutate(ID = ID - 22)


get_functional_data <- function(df, tlen) {
  uid <- unique(df$ID)
  nid <- length(uid)
  
  nk <- (tlen * 100) * (tlen * 100 - 1) / 2
  ji <- df %>%
    group_by(ID) %>%
    summarize(maxtime = max(second)) %>%
    dplyr::select(maxtime) %>%
    unlist() %>%
    unname()
  
  N <- sum(ji)
  smat_w <-
    umat_w <- dmat <- matrix(NA, ncol = nk, nrow = N) # matrix of NAs
  id_vec <- time_g_vec <- rep(NA, N)
  
  inx_samp <- lapply(1:(100 * tlen - 1), function(x) {
    ret <- data.frame("u" = x, "s" = (x + 1):(100 * tlen))
    ret
  })
  inx_samp <- bind_rows(inx_samp)
  inx <- 1
  
  for (i in seq_along(uid)) {
    df_i <- filter(df, ID == i)
    j_i    <- unique(df_i$second)
    for (j in seq_along(j_i)) {
      Y_ij <- df_i$signal_rwrist[df_i$second == j_i[j]]
      umat_w[inx, ] <- Y_ij[inx_samp$u]
      smat_w[inx, ] <- Y_ij[inx_samp$s]
      dmat[inx, ] <- inx_samp$s - inx_samp$u
      time_g_vec[inx] <- ji[j]
      id_vec[inx] <- df_i$ID[1]
      inx <- inx + 1
    }
  }
  
  df_fit <- data.frame(
    "ID" = id_vec,
    "umat" = I(umat_w),
    "smat" = I(smat_w),
    "dmat" = I(dmat),
    "lmat" = I(matrix(1 / nk, ncol = nk, nrow = N))
  )
  rm(list = c("smat_w", "umat_w", "dmat"))
  gc()
  
  
  df_fit <- df_fit %>%
    group_by(ID) %>%
    dplyr::mutate(J = 1:n()) %>%
    ungroup()
  
  df_fit
}

df_fit_zju_s1 <- get_functional_data(df_all_zju_s1, tlen = 1)

df_fit_zju_s1 <- get_functional_data(df_all_zju_s1, tlen = 1)

data_split <- split(df_fit_zju_s1, f = df_fit_zju_s1$ID)

# function to sample percentage
samp <- function(pct, n, ind) {
  set.seed(ind)
  sample(n, floor(pct * n), replace = F)
}
ids <- unique(df_fit_zju_s1$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

getrows <- function(data, rows) {
  subset(data, J %in% rows)
}
getrows_test <- function(data, rows) {
  subset(data, J %notin% rows)
}

data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


fit_functional_model <- function(train, test) {
  uid        <- unique(train$ID)
  nid        <- length(uid)
  
  lpmat <-
    matrix(NA, ncol = nid, nrow = nrow(test)) # to store predictions
  for (i in 1:nid) {
    train$Y <- as.numeric(train$ID == uid[i])
    fit_i <-
      gam(
        Y ~ te(umat, smat, dmat, by = lmat),
        method = "REML",
        data = train,
        family = quasibinomial()
      )
    lp_i <- predict.gam(fit_i, newdata = test, type = "link")
    lpmat[, i] <- unname(lp_i)
    rm(list = c("fit_i", "lp_i"))
    train$Y <- NULL
  }
  lpmat
}

# note this probably needs to be run on a computing cluster
preds <- fit_functional_model(data_train, data_test)


preds %<>%
  data.frame() %>%
  janitor::clean_names() %>%
  expit()

row_sums <- rowSums(preds, na.rm = TRUE)

# normalize and add "true subject column"

preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = data_test$ID)

if (!dir.exists("predictions")) {
  dir.create(here::here("predictions"))
}

saveRDS(preds, here::here("predictions/zjus1_func_predictions.rds"))

# train and test on session 2
df_all_zju_s2 <- read_csv(here::here("data/df_all_zju.csv"),
                          col_types = cols(...1 = col_skip())) %>%
  filter(session == "session_2") %>%
  mutate(ID = ID - 22)

df_fit_zju_s2 <- get_functional_data(df_all_zju_s2, tlen = 1)

data_split <- split(df_fit_zju_s2, f = df_fit_zju_s2$ID)

# function to sample percentage

ids <- unique(df_fit_zju_s2$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)


data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


preds <- fit_functional_model(data_train, data_test)

preds %<>%
  data.frame() %>%
  janitor::clean_names() %>%
  expit()

row_sums <- rowSums(preds)

# normalize and add "true subject column"

preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = data_test$ID)

saveRDS(preds, here::here("predictions/zjus2_func_predictions.rds"))


# train on session 1, test on session 2
df_all_zju <- read_csv(here::here("data/df_all_zju.csv"),
                       col_types = cols(...1 = col_skip())) %>%
  filter(session != "session_0") %>%
  mutate(ID = ID - 22)

df_all_zju_s1 <- read_csv(here::here("data/df_all_zju.csv"),
                          col_types = cols(...1 = col_skip())) %>%
  filter(session == "session_1") %>%
  mutate(ID = ID - 22)

df_all_zju_s2 <- read_csv(here::here("data/df_all_zju.csv"),
                          col_types = cols(...1 = col_skip())) %>%
  filter(session == "session_2") %>%
  mutate(ID = ID - 22)

df_fit_zju_s1 <- get_functional_data(df_all_zju_s1, tlen = 1)
df_fit_zju_s2 <- get_functional_data(df_all_zju_s2, tlen = 1)

df_train <- df_fit_zju_s1
df_test <- df_fit_zju_s2


preds <- fit_functional_model(df_train, df_test)
preds %<>%
  data.frame() %>%
  janitor::clean_names() %>%
  expit()

row_sums <- rowSums(preds)

# normalize and add "true subject column"

preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = df_test$ID)

saveRDS(preds,
        here::here("predictions/zjus1s2_func_predictions.rds"))

# train/test on combined s1 and s2
df_all_zju <- read_csv(here::here("data/df_all_zju.csv"),
                       col_types = cols(...1 = col_skip())) %>%
  filter(session != "session_0") %>%
  mutate(ID = ID - 22)

df_all_zju <-
  df_all_zju %>%
  group_by(ID) %>%
  mutate(s_allrec = row_number(),
         second2 = ceiling(s_allrec / 100)) %>%
  ungroup() %>%
  dplyr::select(-c(second, s_allrec)) %>%
  rename(second = second2)

df_fit_zju <- get_functional_data(df_all_zju, tlen = 1)

data_split <- split(df_fit_zju, f = df_fit_zju$ID)

# function to sample percentage
ids <- unique(df_fit_zju$ID)
# number of rows for each individual
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)


data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


preds <- fit_functional_model(data_train, data_test)

preds %<>%
  data.frame() %>%
  janitor::clean_names() %>%
  expit()

row_sums <- rowSums(preds)

# normalize and add "true subject column"

preds %<>%
  bind_cols(sum = row_sums) %>%
  rowwise() %>%
  mutate(across(x1:x153, ~ .x / sum)) %>%
  dplyr::select(-sum) %>%
  ungroup() %>%
  bind_cols(true_subject = data_test$ID)

saveRDS(preds,
        here::here("predictions/zjuall_func_predictions.rds"))
rm(list = ls())
# compare results 
library(readr)
library(tidyverse)
library(purrr)
library(magrittr)
library(kableExtra)


predictions_filenames <- list.files(here::here("predictions"))
predictions_filenames <- predictions_filenames[grep("predictions", predictions_filenames)]

get_summary <- function(filename){
  tmp <- readRDS(here::here(paste0("predictions/", filename)))
  data <- sub("_.*", "", sub("predictions.rds.*", "", filename))
  method <- str_remove_all(str_extract(filename, "_(.+)_"), "_")
  tmp %<>%
    group_by(true_subject) %>%
    mutate(
      sec = row_number()) %>%
    pivot_longer(cols = -c("true_subject", "sec")) %>%
    mutate(
      model = as.numeric(sub(".*x", "", name))) %>%
    rename(pred = value) %>%
    ungroup() %>%
    group_by(true_subject, model) %>%
    summarize(
      mean_pred = mean(pred, na.rm = TRUE)) %>%
    group_by(true_subject) %>%
    mutate(
      rank = rank(-mean_pred)
    ) %>% 
    filter(model==true_subject) %>%
    mutate(
      rank1 = ifelse(rank == 1, 1, 0),
      rank5 = ifelse(rank <= 5, 1, 0)
    )
  return(data.frame(data = data,
                    method = method,
                    accr1 = sum(tmp$rank1)/nrow(tmp), 
                    accr5 = sum(tmp$rank5)/nrow(tmp),
                    r1 = sum(tmp$rank1), 
                    r5 = sum(tmp$rank5),
                    total = nrow(tmp)))
}

result_summary <-
  map(.x = predictions_filenames,
      .f = get_summary) %>%
  list_rbind()

# make nice table for paper 
result_summary %>%
  filter(data != "zjus2" & data!= "zjuall") %>%
  mutate(across(3:4, ~round(.x, 2))) %>%
  kbl(caption="Summary of Accuracy",
      format="latex",
      col.names = c("Data and Task", "Model", "R1 Accuracy", "R5 Accuracy", 
                    "Rank 1 Correct", "Rank 5 Correct", "Total"),
      align="l") %>%
  kable_minimal(full_width = F,  html_font = "Source Sans Pro")

# compare results averaging over various numbers of seconds in testing data 

get_accuracies <- function(predictions, seconds){
  tmp <- 
    predictions %>%
    group_by(true_subject) %>% 
    mutate(
      sec = floor(row_number()/seconds)) %>% 
    pivot_longer(cols = -c("true_subject", "sec")) %>%
    mutate(
      model = as.numeric(sub(".*x", "", name))) %>%
    rename(pred = value) %>%
    ungroup() %>%
    group_by(true_subject, model, sec) %>%
    summarize(
      mean_pred = mean(pred, na.rm = TRUE)) %>%
    group_by(true_subject, sec) %>%
    mutate(
      rank = rank(-mean_pred)
    ) %>% 
    filter(model==true_subject) %>%
    mutate(
      rank1 = ifelse(rank == 1, 1, 0),
      rank5 = ifelse(rank <= 5, 1, 0)
    )
  n = nrow(tmp)
  r1 = sum(tmp$rank1)/n
  r5 = sum(tmp$rank5)/n
  return(data.frame(metric = c("rank1", "rank5"),
                    value = c(r1, r5),
                    s = c(seconds, seconds)))
  
}


secs <- c(1,2,5,10,15,20,25,30,35,50,60,70,80,90,100)
stats_ml <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/IU_ml_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "ml",
    data = "IU"
  )

stats_logistic <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/IU_logistic_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "logistic",
    data = "IU"
  )

stats_func <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/IU_functional_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "functional",
    data = "IU"
  )

accs_IU <-
  bind_rows(stats_logistic, stats_ml, stats_func) 

secs <- c(1,2,5,10,15,20,25,30,35,50,60,70,80,90,100)
stats_ml <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1_ml_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "ml",
    data = "zjus1"
  )

stats_logistic <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1_logistic_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "logistic",
    data = "zjus1"
  )

stats_func <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1_functional_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "functional",
    data = "zjus1"
  )

accs_zjus1 <-
  bind_rows(stats_logistic, stats_ml, stats_func) 


secs <- c(1,2,5,10,15,20,25,30,35,50,60,70,80,90,100)
stats_ml <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1s2_ml_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "ml",
    data = "zjus1s2"
  )

stats_logistic <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1s2_logistic_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "logistic",
    data = "zjus1s2"
  )

stats_func <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1s2_functional_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "functional",
    data = "zjus1s2"
  )

accs_zjus1s2 <-
  bind_rows(stats_logistic, stats_ml, stats_func) 


all_accs <-
  bind_rows(accs_IU, accs_zjus1, accs_zjus2)

all_accs %>%
  mutate(
    metric = ifelse(metric == "rank1", "Rank-1 Accuracy", "Rank-5 Accuracy")
  ) %>%
  ggplot(aes(x = s, y = value, col = model, linetype = model))+
  geom_line(linewidth=.8, alpha=.9)+
  theme_linedraw()+
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.position =  "bottom",
        legend.text = element_text(size = 12),
        aspect.ratio = 1)+
  facet_grid(metric~data)+
  labs(x = "Number of Seconds", y = "Accuracy", title = "",
       subtitle = "")+
  scale_color_brewer(palette = "Dark2", labels = c("Functional", "Logistic", "ML"), name = "")+
  scale_linetype_manual(name = "", values = c("solid", "dashed", "longdash"), labels = c("Functional", "Logistic", "ML"))

all_accs %>%
  mutate(
    metric = ifelse(metric == "rank1", "Rank-1 Accuracy", "Rank-5 Accuracy")
  ) %>%
  ggplot(aes(x = s, y = value, col = model, shape = model))+
  geom_line(linewidth=.8, alpha=.9)+
  geom_point()+
  theme_linedraw()+
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.position =  "bottom",
        legend.text = element_text(size = 12),
        aspect.ratio = 1)+
  facet_grid(metric~data)+
  labs(x = "Number of Seconds", y = "Accuracy", title = "",
       subtitle = "")+
  scale_color_brewer(palette = "Dark2", labels = c("Functional", "Logistic", "ML"), name = "")+
  scale_shape(name = "", labels = c("Functional", "Logistic", "ML"))



# generate all tables and figures
library(readr)
library(ggplot2)
library(tidyverse)
library(tidymodels)
expit <- function(x) 1/(1+exp(-x))
library(kableExtra)

### figure 1:  problem illustration
df_all_zju <- read_csv("df_all_zju.csv",
                       col_types = cols(...1 = col_skip())) %>%
  filter(session == "session_1") %>%
  mutate(ID = ID - 22)

# theme_set(theme_classic())

set.seed(123)
subs <- sample(unique(df_all_zju$ID), 4, replace = F)
subs <- sort(subs)
supp.labs <-
  c("Subject 14", "Subject 43", "Subject 50", "Subject 118")
names(supp.labs) <- subs
g3 <- df_all_zju %>% filter(ID %in% subs &
                              second <= 3) %>% dplyr::select(signal_rwrist, time, ID) %>%
  ggplot(aes(
    x = time / 100,
    y = signal_rwrist,
    color = as.factor(ID)
  )) + facet_wrap(. ~ ID,
                  nrow = 4,
                  labeller =
                    labeller(ID = supp.labs)) +
  geom_line(linewidth = 1.1) + scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")) +
  theme(legend.position = "none") + scale_y_continuous(limits = c(0, 3),
                                                       breaks = seq(0, 3, 1)) +
  labs(x = "Time (seconds)", y = "Acceleration (g)") + theme_linedraw() + theme(strip.text = element_text(size = 15),
                                                                                axis.text = element_text(size = 12),
                                                                                axis.title = element_text(size = 15),
                                                                                legend.position =  "none")



set.seed(1234)
s2 <- sample(unique(df_all_zju$ID), 2, replace  = F)
subs <- sort(c(118, 43, 28, 80))
supp.labs <- c("Subject ?", "Subject ?", "Subject ?", "Subject ?")
names(supp.labs) <- subs
g4 <- df_all_zju %>% filter(ID %in% c(118, 43, 28, 80) &
                              second <= 6 &
                              second > 3) %>% dplyr::select(signal_rwrist, time, ID) %>%
  ggplot(aes(x = time / 100, y = signal_rwrist, col = as.factor(ID))) + facet_wrap(. ~ ID,
                                                                                   nrow = 4,
                                                                                   labeller =
                                                                                     labeller(ID = supp.labs)) +
  geom_line(linewidth = 1.1) + scale_color_manual(values = c("#66A61E" , "#E6AB02", "#A6761D", "#666666")) +
  theme_linedraw()+ scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  labs(x ="", y = "") +  theme(strip.text = element_text(size = 15),
                               axis.text = element_text(size = 12),
                               axis.title = element_text(size = 15),
                               legend.position =  "none")
g4
ggpubr::ggarrange(g3, g4, common.legend = T, legend = "none")


## figure 2

grid_data_lw_IU <- readRDS("~/Documents/ml_walking_fingerprint/grid_data_lw_IU.rds")
df_all_IU <- read_csv("/Users/lilykoffman/Documents/ml_walking_fingerprint/df_all_IU.csv", 
                      col_types = cols(...1 = col_skip())) %>%
  rename(ID = ID2,
         signal_lwrist = signal_lw,
         signal_lhip = signal_lh,
         signal_lankle = signal_la,
         signal_rankle = signal_ra)

onesig_lag <- function(data, id, s, lag){
  data %>% 
    filter(ID == id) %>%
    filter(second == s) %>%
    mutate(lagsig = lag(signal_lwrist, lag),
           lag = lag)
}

lags <- rep(c(1, 15, 99), 2)
seconds <- c(rep(1, 3), rep(2, 3))

df <-
  map2_dfr(data = df_all_IU, id = 19, 
           .f = onesig_lag,
           .x = seconds,
           .y = lags)
x.labs <- c("u = 1", "u = 15", "u  = 99")
names(x.labs) <- c(1, 15, 99)

y.labs <- c("j = 1", "j = 2")
names(y.labs) <- c(1, 2)

df %>%
  ggplot(aes(x = lagsig, y = signal_lwrist))+
  scale_x_continuous(limits=c(0.75, 1.25), breaks = seq(0.75, 1.25, .25))+
  scale_y_continuous(limits=c(0.75, 1.25), breaks = seq(0.75, 1.25, .25))+
  geom_point()+
  facet_grid(second ~ lag,
             labeller = labeller(second = y.labs, lag = x.labs))+
  theme_linedraw()+
  theme(strip.text = element_text(face = "italic", size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))+
  labs(x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'), title = "")

# figure 3 
df <- 
  df %>% mutate(
    max_x = 1,
    max_y = 1, 
    n = 1
  )
gc_dat <- 
  df %>% 
  group_by(second, lag) %>%
  mutate(cut_x = cut(lagsig, breaks = seq(.75, 1.5, .25), include.lowest = T),
         cut_y = cut(signal_lwrist, breaks = seq(.75, 1.5, .25), include.lowest = T)) %>%
  count(cut_x, cut_y, .drop=FALSE) %>%
  drop_na() %>% 
  mutate(
    cut_x = as.character(cut_x),
    cut_y = as.character(cut_y),
    min_x = str_remove(sub(",.*", "", cut_x), "\\[|\\("),
    max_x = str_remove(sub(".*,", "", cut_x), "\\]|\\)"),
    min_y = str_remove(sub(",.*", "", cut_y), "\\[|\\("),
    max_y = str_remove(sub(".*,", "", cut_y), "\\]|\\)")
  ) %>%
  mutate(across(min_x:max_y, as.numeric)) %>% 
  rename(lagsig = min_x,
         signal_lwrist = min_y)
library(viridis)

p1 <- 
  ggplot(data = gc_dat, aes(xmin  = lagsig, xmax = max_x,
                            ymin = signal_lwrist, ymax = max_y,
                            fill = n))+
  scale_fill_viridis(name = latex2exp::TeX(r'($X_{ijg}$)'))+
  geom_rect(alpha = 0.9)+
  facet_grid(second ~ lag, 
             labeller = labeller(second = y.labs, lag = x.labs))


p2 <- 
  p1 + 
  geom_point(data = df, aes(x = lagsig, y = signal_lwrist), fill = "black")+
  labs(x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), 
       y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'), title = "")+
  scale_x_continuous(limits = c(0.75, 1.5), breaks=seq(.75, 1.5, .25))+
  scale_y_continuous(limits = c(0.75,1.5), breaks=seq(.75, 1.5, .25))+
  geom_hline(yintercept = c(.75, 1, 1.25, 1.5), col = "black")+
  geom_vline(xintercept = c(.75, 1, 1.25, 1.5), col = "black")+
  theme_linedraw()+
  theme(strip.text = element_text(face = "italic", size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

p2 + 
  geom_label(data = gc_dat, 
             aes(x = lagsig + 0.125, y = signal_lwrist + 0.125, label = n),
             col = "white", alpha = 0.8, size = 5)

## figure 5 



df_all_zju <- read_csv("/Users/lilykoffman/Documents/ml_walking_fingerprint/df_all_zju.csv", 
                       col_types = cols(...1 = col_skip())) %>%
  mutate(ID = ID-22)

grid_data_rw_zju_s1 <- readRDS("~/Documents/ml_walking_fingerprint/grid_data_rw_zju_s1.rds") %>%
  mutate(ID = ID - 22)
grid_data_rw_zju_s2 <- readRDS("~/Documents/ml_walking_fingerprint/grid_data_rw_zju_s2.rds") %>%
  mutate(ID = ID - 22)

get_train_test <- function(pct, data){
  data_split <- split(data, f = data$ID)
  
  samp <- function(pct, n, ind) {
    set.seed(ind)
    sample(n, floor(pct * n), replace = F)
  }
  ids <- unique(data$ID)
  # number of rows for each individual 
  rows <- lapply(data_split, nrow) %>% unlist()
  # get random 75% of seconds for training 
  train_indices <- map2(pct = pct,
                        .x = rows,
                        .y = ids,
                        .f = samp)
  
  getrows <- function(data, rows) {
    data[rows, ]
  }
  getrows_test <- function(data, rows) {
    data[-rows, ]
  }
  
  data_train <-
    map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
  data_test <-
    map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)
  return(list(train = data_train, test = data_test))
}

data_train <- get_train_test(pct = .75, data = grid_data_rw_zju_s1)$train
data_test <- get_train_test(pct = .75, data = grid_data_rw_zju_s1)$test



plot_w_points_zju <- function(data_train, data_test, sub, raw_data){
  nzv_trans <- 
    recipe(ID ~ ., data = data_train) %>% 
    step_nzv(all_predictors())
  
  nzv_estimates <- prep(nzv_trans)
  
  nzv <- colnames(juice(nzv_estimates))
  dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv), -second)
  dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv), -second)
  
  train <- dat_nzv
  train$class <- ifelse(train$ID == sub, 1, 0)
  tmp <- train %>% dplyr::select(-c(ID)) 
  mod <- glm(class ~ ., data = tmp, family=binomial(link="logit"))
  
  
  summary <-
    mod %>% 
    tidy() %>%
    arrange(p.value) %>%
    filter(term != "(Intercept)") 
  
  terms <- summary$term
  # variance covariance matrix 
  v_hat <- vcov(mod)[terms, terms]
  a <- 1/(sqrt(diag(v_hat))) 
  A <- a %*% t(a)
  C <- v_hat*A
  
  q <- mvtnorm::qmvnorm(p = .95, corr = C)$quantile
  
  summary <- 
    summary %>% rowwise() %>% 
    mutate(
      lb_marg = estimate - (1.96*std.error),
      ub_marg = estimate + (1.96*std.error),
      lb_joint = estimate - (q*std.error),
      ub_joint = estimate +  (q*std.error),
      sig_marg = ifelse(between(0, lb_marg, ub_marg), 0, 1),
      sig_joint = ifelse(between(0, lb_joint, ub_joint), 0, 1),
    )
  summary <- 
    summary %>% rowwise() %>% 
    mutate(
      lag = str_sub(term, -3, -2),
      sig = sub('.', '', str_split(term, " ")[[1]][1]),
      lagsig = str_split(term, " ")[[1]][2]
    )
  
  if(nrow(summary %>% filter(sig_joint==1)) == 0){
    return(NULL)
  }
  else{
    train_times <-
      data_train %>%
      dplyr::select(ID, second) %>%
      filter(ID == sub)
    
    
    pts_dat <- 
      raw_data %>% inner_join(., train_times, by = c("ID" = "ID", "second" = "second")) %>%
      group_by(ID, second) %>% mutate(
        sig = signal_rwrist, 
        lag_1 = dplyr::lag(signal_rwrist, 15),
        lag_2 = dplyr::lag(signal_rwrist, 30),
        lag_3 = dplyr::lag(signal_rwrist, 45)
      ) %>% dplyr::select(sig, lag_1, lag_2, lag_3, ID, second) %>% 
      pivot_longer(lag_1:lag_3) %>% filter(!is.na(value)) %>% 
      mutate(lag = case_when(
        name == "lag_1" ~ 15,
        name == "lag_2" ~ 30,
        name== "lag_3" ~ 45,
      ),
      xmin = value, 
      xmax = value,
      ymin = sig,
      ymax = sig) 
    
    labs <- c("u = 15", "u = 30", "u = 45")
    names(labs) <- c(15, 30, 45)
    g <- summary %>% 
      mutate(
        est = ifelse(sig_joint == 1, exp(estimate), NA)) %>% 
      mutate(
        xmin = as.numeric(sub(".*\\(", "", sub(",.*", "", lagsig))),
        xmax = as.numeric(str_sub(sub(".*,", "", lagsig), end = -2)),
        ymin =  as.numeric(sub(".*\\(", "", sub(",.*", "", sig))),
        ymax = as.numeric(str_sub(sub(".*,", "", sig), end = -2))
      ) %>% ggplot(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+geom_rect(aes(fill = est)) +
      facet_wrap(.~lag, labeller = labeller(lag = labs))+scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                                                                              midpoint = 1,na.value="white", name = "Odds Ratio")+
      labs(title = "",
           subtitle = "",
           x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'))+
      scale_x_continuous(limits=c(0,3))+scale_y_continuous(limits=c(0,3))+
      theme_linedraw()+
      theme(strip.text = element_text(face = "italic", size = 15),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15))
    
    g + geom_point(data = pts_dat, aes(x = value, y = sig), col = "black", 
                   alpha = 0.03, size = .4)
  }
}

plot_w_points_zju_unadj <- function(data_train, data_test, sub, raw_data){
  nzv_trans <- 
    recipe(ID ~ ., data = data_train) %>% 
    step_nzv(all_predictors())
  
  nzv_estimates <- prep(nzv_trans)
  
  nzv <- colnames(juice(nzv_estimates))
  dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv), -second)
  dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv), -second)
  
  train <- dat_nzv
  train$class <- ifelse(train$ID == sub, 1, 0)
  tmp <- train %>% dplyr::select(-c(ID)) 
  mod <- glm(class ~ ., data = tmp, family=binomial(link="logit"))
  
  
  summary <-
    mod %>% 
    tidy() %>%
    arrange(p.value) %>%
    filter(term != "(Intercept)") 
  
  terms <- summary$term
  # variance covariance matrix 
  v_hat <- vcov(mod)[terms, terms]
  a <- 1/(sqrt(diag(v_hat))) 
  A <- a %*% t(a)
  C <- v_hat*A
  
  q <- mvtnorm::qmvnorm(p = .95, corr = C)$quantile
  
  summary <- 
    summary %>% rowwise() %>% 
    mutate(
      lb_marg = estimate - (1.96*std.error),
      ub_marg = estimate + (1.96*std.error),
      lb_joint = estimate - (q*std.error),
      ub_joint = estimate +  (q*std.error),
      sig_marg = ifelse(between(0, lb_marg, ub_marg), 0, 1),
      sig_joint = ifelse(between(0, lb_joint, ub_joint), 0, 1),
    )
  summary <- 
    summary %>% rowwise() %>% 
    mutate(
      lag = str_sub(term, -3, -2),
      sig = sub('.', '', str_split(term, " ")[[1]][1]),
      lagsig = str_split(term, " ")[[1]][2]
    )
  
  if(nrow(summary %>% filter(sig_marg==1)) == 0){
    return(NULL)
  }
  else{
    train_times <-
      data_train %>%
      dplyr::select(ID, second) %>%
      filter(ID == sub)
    
    
    pts_dat <- 
      raw_data %>% inner_join(., train_times, by = c("ID" = "ID", "second" = "second")) %>%
      group_by(ID, second) %>% mutate(
        sig = signal_rwrist, 
        lag_1 = dplyr::lag(signal_rwrist, 15),
        lag_2 = dplyr::lag(signal_rwrist, 30),
        lag_3 = dplyr::lag(signal_rwrist, 45)
      ) %>% dplyr::select(sig, lag_1, lag_2, lag_3, ID, second) %>% 
      pivot_longer(lag_1:lag_3) %>% filter(!is.na(value)) %>% 
      mutate(lag = case_when(
        name == "lag_1" ~ 15,
        name == "lag_2" ~ 30,
        name== "lag_3" ~ 45,
      ),
      xmin = value, 
      xmax = value,
      ymin = sig,
      ymax = sig) 
    labs <- c("u = 15", "u = 30", "u = 45")
    names(labs) <- c(15, 30, 45)
    g <- summary %>% 
      mutate(
        est = ifelse(sig_marg == 1, exp(estimate), NA)) %>% 
      mutate(
        xmin = as.numeric(sub(".*\\(", "", sub(",.*", "", lagsig))),
        xmax = as.numeric(str_sub(sub(".*,", "", lagsig), end = -2)),
        ymin =  as.numeric(sub(".*\\(", "", sub(",.*", "", sig))),
        ymax = as.numeric(str_sub(sub(".*,", "", sig), end = -2))
      ) %>% ggplot(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+geom_rect(aes(fill = est)) +
      facet_wrap(.~lag, labeller = labeller(lag = labs))+scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                                                                              midpoint = 1,na.value="white", name = "Odds Ratio")+
      labs(title = "",
           subtitle = "",
           x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'))+
      scale_x_continuous(limits=c(0,3))+scale_y_continuous(limits=c(0,3))+
      theme_linedraw()+
      theme(strip.text = element_text(face = "italic", size = 15),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 15))
    
    g + geom_point(data = pts_dat, aes(x = value, y = sig), col = "black", 
                   alpha = 0.03, size = .4)
  }
}


plot_w_points_zju(data_train, data_test, sub = 143, raw_data = df_all_zju)
plot_w_points_zju_unadj(data_train, data_test, sub = 143, raw_data = df_all_zju)


## get legends to match up better 

raw_data <- df_all_zju
sub <- 143

nzv_trans <- 
  recipe(ID ~ ., data = data_train) %>% 
  step_nzv(all_predictors())

nzv_estimates <- prep(nzv_trans)

nzv <- colnames(juice(nzv_estimates))
dat_nzv <- data_train %>% dplyr::select(ID, all_of(nzv), -second)
dat_nzv_test <- data_test %>% dplyr::select(ID, all_of(nzv), -second)

train <- dat_nzv
train$class <- ifelse(train$ID == sub, 1, 0)
tmp <- train %>% dplyr::select(-c(ID)) 
mod <- glm(class ~ ., data = tmp, family=binomial(link="logit"))


summary <-
  mod %>% 
  tidy() %>%
  arrange(p.value) %>%
  filter(term != "(Intercept)") 

terms <- summary$term
# variance covariance matrix 
v_hat <- vcov(mod)[terms, terms]
a <- 1/(sqrt(diag(v_hat))) 
A <- a %*% t(a)
C <- v_hat*A

q <- mvtnorm::qmvnorm(p = .95, corr = C)$quantile

summary <- 
  summary %>% rowwise() %>% 
  mutate(
    lb_marg = estimate - (1.96*std.error),
    ub_marg = estimate + (1.96*std.error),
    lb_joint = estimate - (q*std.error),
    ub_joint = estimate +  (q*std.error),
    sig_marg = ifelse(between(0, lb_marg, ub_marg), 0, 1),
    sig_joint = ifelse(between(0, lb_joint, ub_joint), 0, 1),
  )
summary <- 
  summary %>% rowwise() %>% 
  mutate(
    lag = str_sub(term, -3, -2),
    sig = sub('.', '', str_split(term, " ")[[1]][1]),
    lagsig = str_split(term, " ")[[1]][2]
  )

train_times <-
  data_train %>%
  dplyr::select(ID, second) %>%
  filter(ID == sub)


pts_dat <- 
  raw_data %>% inner_join(., train_times, by = c("ID" = "ID", "second" = "second")) %>%
  group_by(ID, second) %>% mutate(
    sig = signal_rwrist, 
    lag_1 = dplyr::lag(signal_rwrist, 15),
    lag_2 = dplyr::lag(signal_rwrist, 30),
    lag_3 = dplyr::lag(signal_rwrist, 45)
  ) %>% dplyr::select(sig, lag_1, lag_2, lag_3, ID, second) %>% 
  pivot_longer(lag_1:lag_3) %>% filter(!is.na(value)) %>% 
  mutate(lag = case_when(
    name == "lag_1" ~ 15,
    name == "lag_2" ~ 30,
    name== "lag_3" ~ 45,
  ),
  xmin = value, 
  xmax = value,
  ymin = sig,
  ymax = sig) 

labs <- c("u = 15", "u = 30", "u = 45")
names(labs) <- c(15, 30, 45)
max <- max(exp(summary$estimate))
min <- min(exp(summary$estimate))
g <- summary %>% 
  mutate(
    est = ifelse(sig_joint == 1, exp(estimate), NA)) %>% 
  mutate(
    xmin = as.numeric(sub(".*\\(", "", sub(",.*", "", lagsig))),
    xmax = as.numeric(str_sub(sub(".*,", "", lagsig), end = -2)),
    ymin =  as.numeric(sub(".*\\(", "", sub(",.*", "", sig))),
    ymax = as.numeric(str_sub(sub(".*,", "", sig), end = -2))
  ) %>% ggplot(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+geom_rect(aes(fill = est)) +
  facet_wrap(.~lag, labeller = labeller(lag = labs))+scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                                                                          midpoint = 1,na.value="white", name = "Odds Ratio",
                                                                          limits = c(min, max),
                                                                          breaks=seq(0.5, 2, 0.5))+
  labs(title = "",
       subtitle = "",
       x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'))+
  scale_x_continuous(limits=c(0,3))+scale_y_continuous(limits=c(0,3))+
  theme_linedraw()+
  theme(strip.text = element_text(face = "italic", size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        aspect.ratio = 1,
        legend.position = "bottom")

g + geom_point(data = pts_dat, aes(x = value, y = sig), col = "black", 
               alpha = 0.03, size = .4)

g <- summary %>% 
  mutate(
    est = ifelse(sig_marg == 1, exp(estimate), NA)) %>% 
  mutate(
    xmin = as.numeric(sub(".*\\(", "", sub(",.*", "", lagsig))),
    xmax = as.numeric(str_sub(sub(".*,", "", lagsig), end = -2)),
    ymin =  as.numeric(sub(".*\\(", "", sub(",.*", "", sig))),
    ymax = as.numeric(str_sub(sub(".*,", "", sig), end = -2))
  ) %>% ggplot(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+geom_rect(aes(fill = est)) +
  facet_wrap(.~lag, labeller = labeller(lag = labs))+scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                                                                          midpoint = 1,na.value="white", name = "Odds Ratio",
                                                                          limits = c(min, max),
                                                                          breaks=seq(0.5, 2, 0.5))+
  labs(title = "",
       subtitle = "",
       x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'))+
  scale_x_continuous(limits=c(0,3))+scale_y_continuous(limits=c(0,3))+
  theme_linedraw()+
  theme(strip.text = element_text(face = "italic", size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "bottom",
        aspect.ratio = 1)

g + geom_point(data = pts_dat, aes(x = value, y = sig), col = "black", 
               alpha = 0.03, size = .4)

# fingerprints (figure 6?)

onelag <- function(subject, lag, data){
  df_dens <- 
    data %>% 
    filter(ID==subject) %>% 
    dplyr::select(signal_rwrist, second) %>% 
    group_by(second) %>%
    mutate(
      lag_signal = lag(signal_rwrist, n = lag)) %>% 
    drop_na() %>% 
    mutate(
      lag = lag 
    )
  df_dens$density <- get_density(df_dens$signal_rwrist, df_dens$lag_signal, n=100)
  df_dens
}
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


df_all_zju <- read_csv("/Users/lilykoffman/Documents/ml_walking_fingerprint/df_all_zju.csv", 
                       col_types = cols(...1 = col_skip())) %>%
  mutate(
    ID = ID - 22
  )
df_zju_train <-
  df_all_zju %>% 
  filter(session == "session_1") 

df_zju_test <-
  df_all_zju %>% 
  filter(session == "session_2")
library(viridis)

labs <- c("u = 15", "u = 30", "u = 45")
names(labs) <- c(15, 30, 45)
all_densities_train <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_train, subject=5) %>%
  mutate(session = "Session 1")
all_densities_test <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_test, subject=5)  %>%
  mutate(session = "Session 2")

all_dens <- bind_rows(all_densities_test, all_densities_train)
maxdens <- max(all_dens$density)

s5 <- 
  all_dens %>% 
  ggplot(aes(x = lag_signal, y = signal_rwrist, col = density))+
  geom_point() + 
  scale_color_viridis(name = "Density", limits=c(0, maxdens))+
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  theme_linedraw() + facet_grid(session~lag, labeller = labeller(lag = labs))+ 
  labs(x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), 
       y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'),
       title = "Subject 5")+
  theme(strip.text.x = element_text(size = 15, face = "italic"),
        strip.text.y = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16),
        aspect.ratio = 1,
        legend.position = "bottom")
s5
all_densities_train <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_train, subject=79) %>%
  mutate(session = "Session 1")
all_densities_test <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_test, subject=79)  %>%
  mutate(session = "Session 2")

all_dens <- bind_rows(all_densities_test, all_densities_train)
maxdens <- max(all_dens$density)

s79 <- 
  all_dens %>% 
  ggplot(aes(x = lag_signal, y = signal_rwrist, col = density))+
  geom_point() + 
  scale_color_viridis(name = "Density", limits=c(0, maxdens))+
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  theme_linedraw() + facet_grid(session~lag, labeller = labeller(lag = labs))+ 
  labs(x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), 
       y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'),
       title = "Subject 79")+
  theme(strip.text.x = element_text(size = 15, face = "italic"),
        strip.text.y = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16),
        aspect.ratio = 1,
        legend.position = "bottom")

s79
ggpubr::ggarrange(s5, s79)

all_densities_train <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_train, subject=3) %>%
  mutate(session = "Session 1")
all_densities_test <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_test, subject=3)  %>%
  mutate(session = "Session 2")

all_dens <- bind_rows(all_densities_test, all_densities_train)
maxdens <- max(all_dens$density)

s3 <- 
  all_dens %>% 
  ggplot(aes(x = lag_signal, y = signal_rwrist, col = density))+
  geom_point() + 
  scale_color_viridis(name = "Density", limits=c(0, maxdens))+
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  theme_linedraw() + facet_grid(session~lag, labeller = labeller(lag = labs))+ 
  labs(x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), 
       y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'),
       title = "Subject 3")+
  theme(strip.text.x = element_text(size = 15, face = "italic"),
        strip.text.y = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16),
        aspect.ratio = 1,
        legend.position = "bottom")
all_densities_train <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_train, subject=136) %>%
  mutate(session = "Session 1")
all_densities_test <- map_dfr(.x = seq(15,45,15), .f = onelag, data = df_zju_test, subject=136)  %>%
  mutate(session = "Session 2")

all_dens <- bind_rows(all_densities_test, all_densities_train)
maxdens <- max(all_dens$density)

s136 <- 
  all_dens %>% 
  ggplot(aes(x = lag_signal, y = signal_rwrist, col = density))+
  geom_point() + 
  scale_color_viridis(name = "Density", limits=c(0, maxdens))+
  scale_x_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1), limits = c(0,3),
                     minor_breaks = seq(0, 3, 0.25)) +
  theme_linedraw() + facet_grid(session~lag, labeller = labeller(lag = labs))+ 
  labs(x = latex2exp::TeX(r'(Acceleration at \it{$v(s-u)$} (g))'), 
       y = latex2exp::TeX(r'(Acceleration at \it{$v(s)$} (g))'),
       title = "Subject 136")+
  theme(strip.text.x = element_text(size = 15, face = "italic"),
        strip.text.y = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16),
        aspect.ratio = 1,
        legend.position = "bottom")
ggpubr::ggarrange(s3, s136)

## accuracies figure 
get_accuracies <- function(predictions, seconds){
  tmp <- 
    predictions %>%
    group_by(true_subject) %>% 
    mutate(
      sec = floor(row_number()/seconds)) %>% 
    pivot_longer(cols = -c("true_subject", "sec")) %>%
    mutate(
      model = as.numeric(sub(".*x", "", name))) %>%
    rename(pred = value) %>%
    ungroup() %>%
    group_by(true_subject, model, sec) %>%
    summarize(
      mean_pred = mean(pred, na.rm = TRUE)) %>%
    group_by(true_subject, sec) %>%
    mutate(
      rank = rank(-mean_pred)
    ) %>% 
    filter(model==true_subject) %>%
    mutate(
      rank1 = ifelse(rank == 1, 1, 0),
      rank5 = ifelse(rank <= 5, 1, 0)
    )
  n = nrow(tmp)
  r1 = sum(tmp$rank1)/n
  r5 = sum(tmp$rank5)/n
  return(data.frame(metric = c("rank1", "rank5"),
                    value = c(r1, r5),
                    s = c(seconds, seconds)))
  
}


secs <- c(1,2,5,10,15,20,25,30,35,50,60,70,80,90,100)
stats_ml <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/IU_ml_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "ml",
    data = "IU"
  )

stats_logistic <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/IU_logistic_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "logistic",
    data = "IU"
  )

stats_func <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/IU_func_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "functional",
    data = "IU"
  )

accs_IU <-
  bind_rows(stats_logistic, stats_ml, stats_func) 

secs <- c(1,2,5,10,15,20,25,30,35,50,60,70,80,90,100)
stats_ml <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1_ml_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "ml",
    data = "zjus1"
  )

stats_logistic <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1_logistic_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "logistic",
    data = "zjus1"
  )

stats_func <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1_func_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "functional",
    data = "zjus1"
  )

accs_zjus1 <-
  bind_rows(stats_logistic, stats_ml, stats_func) 


secs <- c(1,2,5,10,15,20,25,30,35,50,60,70,80,90,100)
stats_ml <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1s2_ml_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "ml",
    data = "zjus1s2"
  )

stats_logistic <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1s2_logistic_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "logistic",
    data = "zjus1s2"
  )

stats_func <- 
  map(.x = secs,
      .f = get_accuracies,
      predictions = readRDS(here::here("predictions/zjus1s2_func_predictions.rds"))) %>%
  list_rbind() %>%
  mutate(
    model = "functional",
    data = "zjus1s2"
  )

accs_zjus1s2 <-
  bind_rows(stats_logistic, stats_ml, stats_func) 


all_accs <-
  bind_rows(accs_IU, accs_zjus1, accs_zjus1s2)

all_accs %>%
  mutate(
    metric = ifelse(metric == "rank1", "Rank-1", "Rank-5"),
    data = case_when(
      data == "zjus1" ~ "ZJU S1",
      data == "zjus1s2" ~ "ZJU S1S2",
      TRUE ~ data
    )
  ) %>%
  ggplot(aes(x = s, y = value, col = model, linetype = model))+
  geom_line(linewidth=.9, alpha=1)+
  theme_linedraw()+
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.position =  "bottom",
        legend.text = element_text(size = 12),
        aspect.ratio = 1)+
  facet_grid(metric~data)+
  labs(x = "Number of Seconds", y = "Accuracy", title = "",
       subtitle = "")+
  scale_color_brewer(palette = "Dark2", labels = c("Functional", "Logistic", "ML"), name = "")+
  scale_linetype_manual(name = "", values = c("solid", "dashed", "longdash"), labels = c("Functional", "Logistic", "ML"))

all_accs %>%
  mutate(
    metric = ifelse(metric == "rank1", "Rank-1", "Rank-5")
  ) %>%
  ggplot(aes(x = s, y = value, col = model, shape = model))+
  geom_line(linewidth=.8, alpha=.9)+
  geom_point()+
  theme_linedraw()+
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.position =  "bottom",
        legend.text = element_text(size = 12),
        aspect.ratio = 1)+
  facet_grid(metric~data)+
  labs(x = "Number of Seconds", y = "Accuracy", title = "",
       subtitle = "")+
  scale_color_brewer(palette = "Dark2", labels = c("Functional", "Logistic", "ML"), name = "")+
  scale_shape(name = "", labels = c("Functional", "Logistic", "ML"))

# table 1 


# summary of IU data 
grid_data_lw_IU <- readRDS(here::here("data/grid_data_lw_IU.rds"))

# split into training and testing
# 75% train, 25% test, equal proportions for ea person 
data_split <- split(grid_data_lw_IU, f = grid_data_lw_IU$ID)
# function to sample percentage 
samp <- function(pct, n, ind) {
  set.seed(ind)
  sample(n, floor(pct * n), replace = F)
}
ids <- unique(grid_data_lw_IU$ID)
# number of rows for each individual 
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training 
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)

getrows <- function(data, rows) {
  data[rows, ]
}
getrows_test <- function(data, rows) {
  data[-rows, ]
}

data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


IU_summary <- 
  data_train %>% 
  group_by(ID) %>%
  summarize(n_seconds = n()) %>%
  mutate(
    type = "train"
  ) %>%
  bind_rows(data_test %>%
              group_by(ID) %>%
              summarize(n_seconds = n()) %>%
              mutate(
                type = "test"
              )) %>%
  pivot_wider(names_from = type, values_from = n_seconds) %>%
  ungroup() %>%
  summarize(
    n_test = sum(test),
    n_train = sum(train),
    avg_test = median(test),
    avg_train = median(train),
    iqr_test = IQR(test),
    iqr_train = IQR(train)) %>%
  mutate(data = "IU")

# summary of zju data 
grid_data_rw_zju_s1 <- readRDS(here::here("data/grid_data_rw_zju_s1.rds")) %>%
  mutate(ID = ID - 22)
grid_data_rw_zju_s2 <- readRDS(here::here("data/grid_data_rw_zju_s2.rds")) %>%
  mutate(ID = ID - 22)

data_split <- split(grid_data_rw_zju_s1, f = grid_data_rw_zju_s1$ID)
# function to sample percentage 

ids <- unique(grid_data_rw_zju_s1$ID)
# number of rows for each individual 
rows <- lapply(data_split, nrow) %>% unlist()
# get random 75% of seconds for training 
train_indices <- map2(pct = .75,
                      .x = rows,
                      .y = ids,
                      .f = samp)
data_train <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows)
data_test <-
  map2_dfr(.x = data_split, .y = train_indices, .f = getrows_test)


zju1_summary <- 
  data_test %>% 
  group_by(ID) %>%
  summarize(n_seconds = n()) %>%
  mutate(
    type = "test"
  ) %>%
  bind_rows(data_train %>%
              group_by(ID) %>%
              summarize(n_seconds = n()) %>%
              mutate(
                type = "train"
              )) %>%
  pivot_wider(names_from = type, values_from = n_seconds) %>%
  ungroup() %>%
  summarize(
    n_test = sum(test),
    n_train = sum(train),
    avg_test = median(test),
    avg_train = median(train),
    iqr_test = IQR(test),
    iqr_train = IQR(train)) %>%
  mutate(data = "zjus1")


data_train <- grid_data_rw_zju_s1
data_test <- grid_data_rw_zju_s2



zju12_summary <- 
  data_test %>% 
  group_by(ID) %>%
  summarize(n_seconds = n()) %>%
  mutate(
    type = "test"
  ) %>%
  bind_rows(data_train %>%
              group_by(ID) %>%
              summarize(n_seconds = n()) %>%
              mutate(
                type = "train"
              )) %>%
  pivot_wider(names_from = type, values_from = n_seconds) %>%
  ungroup() %>%
  summarize(
    n_test = sum(test),
    n_train = sum(train),
    avg_test = median(test),
    avg_train = median(train),
    iqr_test = IQR(test),
    iqr_train = IQR(train)) %>%
  mutate(data = "zjus1s2")


table <- 
  bind_rows(IU_summary, zju1_summary, zju12_summary) %>%
  dplyr::select(data,  avg_train, avg_test, iqr_train, iqr_test)

table %>%
  kbl(caption="Summary of Minutes Used for Training and Testing",
      format="latex",
      col.names = c("Data",  "Median Seconds Training", 
                    "Median Seconds Testing", "IQR Training", "IQR Testing"),
      align="l") %>%
  kable_minimal(full_width = F,  html_font = "Source Sans Pro")

# table 2
predictions_filenames <- list.files(here::here("predictions"))
predictions_filenames <- predictions_filenames[grep("predictions", predictions_filenames)]

get_summary <- function(filename){
  tmp <- readRDS(here::here(paste0("predictions/", filename)))
  data <- sub("_.*", "", sub("predictions.rds.*", "", filename))
  method <- str_remove_all(str_extract(filename, "_(.+)_"), "_")
  tmp %<>%
    group_by(true_subject) %>%
    mutate(
      sec = row_number()) %>%
    pivot_longer(cols = -c("true_subject", "sec")) %>%
    mutate(
      model = as.numeric(sub(".*x", "", name))) %>%
    rename(pred = value) %>%
    ungroup() %>%
    group_by(true_subject, model) %>%
    summarize(
      mean_pred = mean(pred, na.rm = TRUE)) %>%
    group_by(true_subject) %>%
    mutate(
      rank = rank(-mean_pred)
    ) %>% 
    filter(model==true_subject) %>%
    mutate(
      rank1 = ifelse(rank == 1, 1, 0),
      rank5 = ifelse(rank <= 5, 1, 0)
    )
  return(data.frame(data = data,
                    method = method,
                    accr1 = sum(tmp$rank1)/nrow(tmp), 
                    accr5 = sum(tmp$rank5)/nrow(tmp),
                    r1 = sum(tmp$rank1), 
                    r5 = sum(tmp$rank5),
                    total = nrow(tmp)))
}

result_summary <-
  map(.x = predictions_filenames,
      .f = get_summary) %>%
  list_rbind()

# make nice table for paper 
result_summary %>%
  filter(data != "zjus2" & data!= "zjuall") %>%
  mutate(across(3:4, ~round(.x, 2))) %>%
  kbl(caption="Summary of Accuracy",
      format="latex",
      col.names = c("Data and Task", "Model", "R1 Accuracy", "R5 Accuracy", 
                    "Rank 1 Correct", "Rank 5 Correct", "Total"),
      align="l") %>%
  kable_minimal(full_width = F,  html_font = "Source Sans Pro")



