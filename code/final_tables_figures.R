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
  