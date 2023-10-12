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
