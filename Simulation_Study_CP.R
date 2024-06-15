############################################################
rm(list=ls())
library(parallel)
library(ggplot2)
library(xgboost)
library(tidymodels)
library(modeltime)
library(tidyverse)
library(timetk)
library(reshape2)
setwd("~/CPTS/")
interactive<-F
############################################################
# beta <- c(-140, -2.73e-04,  1.83e-03,  3.55e-04 ,-3.80e-04 )
# Xmat <- read.csv("./simul.csv",row.names = 1)
# Xmat <- cbind(index=as.character(1:nrow(Xmat)),Xmat)
# df_tmp <- melt(Xmat)
# df_tmp$index <- as.numeric(df_tmp$index)
# ggplot(data=df_tmp,aes(x=index,y=value,group=variable,color=variable)) +
#   geom_point() + geom_line()

# Xmat <- read.csv("./simul.csv",row.names = 1)
# Xmat <- as.matrix(Xmat)

n <- 72
beta <- c(1,0.5,-0.5,1,-1)

RES <- mclapply(1:1000,function(sim_num){
  set.seed(123 + sim_num)
  
  Xmat <- matrix(rnorm(n*(length(beta)-1)),nrow=n)
  dim(Xmat)
  epsvec <- arima.sim(n-1, model = list(order = c(0,1,0)))
  
  y <- cbind(1,Xmat)%*%beta + epsvec
  
  df_tr_val0 <- data.frame(y = y, Xmat)[1:48,]
  lm_result <- lm(y ~ ., data = df_tr_val0)
  
  residuals <- as.vector(lm_result$residuals)
  y_pred <- predict(lm_result,data.frame(Xmat[49:n,]))
  y_true <- y[49:n]
  
  date_vec <- seq(as.Date("2020-01-01"), as.Date("2023-12-01"), by = "month")
  date_vec_future <- seq(as.Date("2024-01-01"), as.Date("2025-12-01"), by = "month")
  
  df_tr_val <- data.frame(Date = date_vec, Values = residuals)
  
  splits <- initial_time_split(df_tr_val, prop = 0.9)
  
  suppressMessages(model_fit_arima_no_boost <- arima_reg() %>%
                     set_engine(engine = "auto_arima") %>%
                     fit(Values ~ Date, data = training(splits)))
  
  suppressMessages(model_fit_arima_boosted <- arima_boost(min_n = 2, learn_rate = 0.1) %>%
                     set_engine(engine = "auto_arima_xgboost") %>%
                     fit(Values ~ Date + as.numeric(Date) + 
                           factor(month(Date, label = TRUE),ordered = FALSE),
                         data = training(splits)))
  
  suppressMessages(model_fit_ets <- exp_smoothing() %>%
                     set_engine(engine = "ets") %>%
                     fit(Values ~ Date, data = training(splits)))
  
  suppressMessages(model_fit_prophet <- prophet_reg() %>%
                     set_engine(engine = "prophet") %>%
                     fit(Values ~ Date, data = training(splits)))
  
  models_tbl <- modeltime_table(
    model_fit_arima_no_boost,
    model_fit_arima_boosted,
    model_fit_ets,
    model_fit_prophet
  )
  
  calibration_tbl <- models_tbl %>%
    modeltime_calibrate(new_data = testing(splits))
  
  
  new_data <- data.frame(Date = date_vec_future)
  
  
  # forecast_tbl <- calibration_tbl %>%
  #   modeltime_forecast(new_data = testing(splits), actual_data = df_tr_val,
  #                      conf_interval = 0.90, conf_method = "conformal_default",
  #      keep_data = TRUE)
  #View(forecast_tbl)
  
  refit_tbl <- calibration_tbl %>%
    modeltime_refit(df_tr_val)
  
  
  
  forecast_tbl <- refit_tbl %>%
    modeltime_forecast(
      new_data      = new_data,
      actual_data   = df_tr_val,
      conf_interval = 0.90,
      conf_method   = "conformal_default", 
      keep_data     = TRUE
    )
  
  
  #View(forecast_tbl)
  
  forecast_tbl <- forecast_tbl %>%
    filter(.model_id != "NA")
  
  patient <- y_true
  fv <- y_pred
  
  res <- lapply(1:4,function(ind){
    ind_vec <- 1:24 + (ind-1)*24
    res_tmp <- (forecast_tbl$.conf_lo[ind_vec] <= y_true-y_pred)&
      (y_true-y_pred <= forecast_tbl$.conf_hi[ind_vec])
    return(res_tmp)
  })
  res <- do.call("cbind",res)
  res <- as.data.frame(res)
  res <- cbind(date_vec_future,res)
  colnames(res) <- c("Date",paste0("m",1:4))
  res$sim_num <- sim_num
  return(res)
},mc.cores = 8)

RES <- do.call("rbind",RES)

RES_trans <- RES %>%
  group_by(Date) %>%
  summarise(Arima = mean(m1), Arima_Boosted = mean(m2),
            ETS = mean(m3), Prophet = mean(m4))
df_tmp <- melt(RES_trans,id.vars="Date")
dim(df_tmp)
ggplot(data=df_tmp,aes(x=Date,y=value,group=variable,color=variable)) +
  geom_point() + geom_line() + theme_bw()


