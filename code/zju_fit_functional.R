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
