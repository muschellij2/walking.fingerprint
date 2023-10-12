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
