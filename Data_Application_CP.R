##############################
rm(list=ls())
#library
library(ggplot2)
library(xgboost)
library(tidymodels)
library(modeltime)
library(tidyverse)
library(timetk)
library(reshape2)
library(forecast)
setwd("~/CPTS/")
#############################
#data

#Prescription dat + the number of each disease
merged_dat<-read.csv("merged_dat.csv")
dim(merged_dat)
str(merged_dat)
colnames(merged_dat)
View(merged_dat)
rownames(merged_dat)<-merged_dat[,1]
merged_dat<-merged_dat[,-1]
##############################
interactive<-F

#coverage rate 
conf_rate<-c() #prophet model 
conf_rate_2<-c() #arima without CP process 
dates<-seq(as.Date("2020-04-01"),as.Date("2023-04-01"), by="1 month")
future_dates <- seq(as.Date("2022-07-01"), as.Date("2023-04-01"), by = "month")
fitted_dat<-data.frame(date=future_dates)
i
colnames(merged_dat)
for (i in 1:15){
  lm_result <- lm(merged_dat[,i]~atc1+atc2+atc3+atc4,data = merged_dat) 
  values<-lm_result$residuals
  df<-data.frame(Date=dates,Values=values)
  #data split 9:1
  splits <- initial_time_split(df, prop = 0.9)
  
  
  # Model 4: prophet ----
  model_fit_prophet <- prophet_reg() %>%
    set_engine(engine = "prophet") %>%
    fit(Values ~ Date, data = training(splits))
  
  #models table 
  models_tbl <- modeltime_table(
    model_fit_prophet)
  
  
  
  #calibrate
  calibration_tbl <- models_tbl %>%
    modeltime_calibrate(new_data = testing(splits))
  
  # # calibration_tbl %>%
  #   modeltime_accuracy() %>%
  #   table_modeltime_accuracy(
  #     .interactive = interactive
  #   )
  

  new_data <- data.frame(
    Date = future_dates
  )
  
  forecast_tbl <- calibration_tbl %>%
    modeltime_forecast(
      new_data      = new_data,  
      actual_data   = df,
      conf_interval = 0.90,
      conf_method   = "conformal_default",
      keep_data     = TRUE
    )
  
  ###########################################3/25
  
  refit_tbl <- calibration_tbl %>%
    modeltime_refit(df)
  
  
  forecast_future_tbl <- refit_tbl %>%
    modeltime_forecast(
      new_data      = new_data,
      actual_data   = df,
      conf_interval = 0.90,
      conf_method   = "conformal_default", 
      keep_data     = TRUE
    )
  future_tbl <- forecast_future_tbl %>%
    filter(.model_id == 1)
  
  fv<-lm_result$fitted.values
  low_cp<-fv[28:37]+future_tbl$.conf_lo
  hi_cp<-fv[28:37]+future_tbl$.conf_hi
  patient<-merged_dat[c(28:37),i]
  future_tbl<-cbind(future_tbl,low_cp,hi_cp,patient)
  
  #ARIMA
  #library(forecast)
  ts <- ts(df$Values[1:27],start = c(2020, 4),end=c(2022,6),frequency = 12 )
  arima_model <- auto.arima(ts)
  fore_arima<-forecast(arima_model,h=10)
  
  low_Arima<-fv[28:37]+fore_arima$lower[,2]
  hi_Arima<-fv[28:37]+fore_arima$upper[,2]
  future_tbl<-cbind(future_tbl,low_Arima,hi_Arima)
  
  #coverage rate
  future_tbl <- future_tbl %>%
    mutate(between_conf = ifelse(patient > low_cp & patient < hi_cp, TRUE, FALSE))
  
  future_tbl<- future_tbl %>%
    mutate(between_conf_arima = ifelse(patient > low_Arima & patient < hi_Arima, TRUE, FALSE))
  forecast_results <- sum(future_tbl$between_conf)
  temp<-{sum(forecast_results)}/10
  conf_rate<-append(conf_rate,temp)
  forecast_results_2 <- sum(future_tbl$between_conf_arima)
  temp<-{sum(forecast_results_2)}/10
  conf_rate_2<-append(conf_rate_2,temp)
  
  temp<-forecast_future_tbl[28:37,ncol(forecast_future_tbl)] + fv[28:37]
  fitted_dat<-cbind(fitted_dat,temp)
  
  combined_fitted_data <- data.frame(date=dates_23,low=c(rep(NA, length(dates_23) - 10),low_cp),hi=c(rep(NA, length(dates_23) - 10),hi_cp),actual=merged_dat[,i],
                                     low_Arima=c(rep(NA, length(dates_23) - 10),low_Arima),hi_Arima=c(rep(NA, length(dates_23) - 10),hi_Arima))
  directory <- "//Users//seheekim//R_2023//sehee_project_result"
  file_path <- file.path(directory, paste("result", colnames(merged_dat)[i], ".csv", sep = "_"))
  write.csv(future_tbl, file = file_path)
  
  
  
  #library(ggplot2)
  

  fill_data <- combined_fitted_data[28:37, ]
  fill_data$date <- as.Date(fill_data$date)  

  fill_poly <- data.frame(
    date = c(fill_data$date, rev(fill_data$date)),
    value = c(fill_data$low, rev(fill_data$hi))
  )
  
  
  # 그래프를 그리는 함수
  plot_graph <- function(data, title) {
    ggplot(data, aes(x = date)) +
      geom_line(aes(y = low, color = "Low")) +
      geom_line(aes(y = hi, color = "High")) +
      geom_line(aes(y = actual, color = "Actual")) +
      geom_line(aes(y = low_Arima, color = "Low ARIMA"), linetype = "dashed") +  # ARIMA 예측값의 하한선
      geom_line(aes(y = hi_Arima, color = "High ARIMA"), linetype = "dashed") +  # ARIMA 예측값의 상한선
      labs(title = title, x = "Date", y = "Value") +
      scale_color_manual(values = c("Low" = "blue", "High" = "red", "Actual" = "black", "Low ARIMA" = "blue", "High ARIMA" = "red")) +
      theme_minimal()
  }
  

  save_plot <- function(graph, file_path, width = 10, height = 5, dpi = 300) {
    ggsave(filename = file_path, plot = graph, width = width, height = height, dpi = dpi)
  }
  

  graph <- plot_graph(combined_fitted_data, paste("Combined Plot", colnames(merged_dat)[i]))
  
  

  directory <- "//Users//seheekim//R_2023//sehee_project_result"
  file_path <- file.path(directory, paste("plot", colnames(merged_dat)[i], ".png", sep = "_"))
  save_plot(graph, file_path)
}


coverage_rate<-data.frame(disease=colnames(merged_dat[,1:15]),CP=conf_rate,ARIMA=conf_rate_2)
length(conf_rate)
length(conf_rate_2)

# Melt the data frame for easy plotting
coverage_rate_long <- tidyr::pivot_longer(coverage_rate, cols = c("CP", "ARIMA"), names_to = "Method", values_to = "Coverage_Rate")

# Plot
ggplot(coverage_rate_long, aes(x = disease, y = Coverage_Rate, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Disease", y = "Coverage Rate", fill = "Method") +
  ggtitle("Coverage Rate Comparison between CP and ARIMA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability

calibration_tbl %>%
  modeltime_accuracy() %>%
  table_modeltime_accuracy(
    .interactive = interactive
  )
