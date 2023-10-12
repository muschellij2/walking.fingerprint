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


