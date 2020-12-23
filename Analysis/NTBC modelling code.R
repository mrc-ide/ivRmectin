
rm(list=ls())

setwd("~/OneDrive - LSTM/ivRmectin")

devtools::test()
Yes
sessionInfo()

R.version



install.packages("devtools") # (if needed)
devtools::install_gitub("mrc-ide/dde", upgrade = FALSE)
                                                                                                                                                                                       install.packages("ggplot2")
install.packages("deSolve")
library(deSolve)
install.packages("RColorBrewer")
install.packages("plyr")
install.packages("dplyr")
install.packages("tidyverse")

odin::odin_package(".")

as.data.frame(installed.packages()[,c(1,3,4)])

install.packages("dde")
library(dde)

use_dde
rm(list = c("effective_mortality", "ivm_fun", "mda_fun", "new_params"))
Yes
# Loading the ivRmectin package
devtools::load_all()
devtools::load_all()

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(tidyverse)

install.packages("drat") # -- if you don't have drat installed
drat:::add("mrc-ide")
install.packages("odin")

devtools::install_github("mrc-ide/odin", upgrade = TRUE)


devtools::install_github("mrc-ide/dde", force = TRUE)
devtools::document()


# Create a vector of age categories for the model
init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 25

#low - 2, moderate - 15, high - 120

# Provide the length of time (in days) that you want to run the model for
time_period <- 3*365 # run model for 10 years

getwd()

set
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

rm(list = c("effective_mortality", "ivm_fun", "mda_fun", "new_params"))


# Load ivermectin hazardd



IVM_3_300_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_IVM_3_300_vitro.RDS")
IVM_3_600_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_IVM_3_600_vitro.RDS")
IVM_1_400_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_IVM_1_400_vitro.RDS")
IVM_3_300_vivo = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_IVM_3_300_vivo.RDS")
IVM_3_600_vivo = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_IVM_3_600_vivo.RDS")
IVM_1_400_vivo = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_IVM_1_400_vivo.RDS")

NTBC_1_100_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_1_100_vitro.RDS")
NTBC_3_100_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_3_100_vitro.RDS")
NTBC_1_300_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_1_300_vitro.RDS")
NTBC_3_300_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_3_300_vitro.RDS")
NTBC_1_600_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_1_600_vitro.RDS")
NTBC_3_600_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_3_600_vitro.RDS")
NTBC_1_1mg_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_1_1mg_vitro.RDS")
NTBC_3_1mg_vitro = readRDS("/Users/anna.trett/Dropbox/My Mac (AdminMac16’s MacBook Pro)/Downloads/IVMHR/01-Data/HRprofiles_NTBC_3_1mg_vitro.RDS")


list_all = list (IVM_3_300_vitro, IVM_3_600_vitro, IVM_1_400_vitro, IVM_3_300_vivo, IVM_3_600_vivo, IVM_1_400_vivo,NTBC_1_100_vitro,
                 NTBC_3_100_vitro,NTBC_1_300_vitro,NTBC_3_300_vitro,NTBC_1_600_vitro,NTBC_3_600_vitro,NTBC_1_1mg_vitro)

# This function runs the odin model with all the parameters that you previously generated with
# the ivRmectimn::create_r_model call
runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}


#loop to obtain optimal timing of seasonal MDA
#warning - takes days to run
##Function - Optimisation


opt_date_Func = function (HR_data) {

  try <- ldply(IVM_3_300_vitro, data.frame)
  med_clin = ddply(.data = try, .variables = "TIME", summarise, "HR" = median(HR))

  data_frame_list  = list()

  sum_inc_list= list()


  vector=vector(mode="numeric",length=3)
  list=vector(mode="list")

  for (a in 1:306) {
    vector[1] <- a
    vector[2] <- vector[1] + 30
    vector[3] <- vector[2] + 30
    list[[a]] <- vector

  }

  for (j in 1:306) {
    init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
    init_EIR <- 25

    time_period <- 365*5# run model for 5 years
    lower_time <- 365*3
    vector<-list[[1]]
    print(j)
    administration = unlist(vector)

    administration[1] = administration[1] + lower_time
    administration[2] = administration[2] + lower_time
    administration[3] = administration[3] + lower_time
    ivm_parms0 = ivm_fun(IVM_start_times = administration,  # time endectocide delivery occurs
                         time_period = time_period, # time period for the model to run over
                         hazard_profile = med_clin[,2], # dummy hazard profile - must be vector
                         ivm_coverage = 0.8, # proportion of population receiving the endectocide
                         ivm_min_age = 5, # youngest age group receiving endectocide
                         ivm_max_age = 90) # oldest age group receiving endectocide

    wh <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                     num_int = 2,
                                     het_brackets = 5,
                                     age = init_age,
                                     init_EIR = 25,
                                     country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                     admin2 = "Fatick",
                                     ttt = ivm_parms0$ttt,
                                     eff_len = ivm_parms0$eff_len,
                                     haz = ivm_parms0$haz,
                                     ivm_cov_par = ivm_parms0$ivm_cov_par,
                                     ivm_min_age = ivm_parms0$ivm_min_age,
                                     ivm_max_age = ivm_parms0$ivm_max_age,
                                     IVRM_start = ivm_parms0$IVRM_start)

    res <- runfun(wh)
    res = as.data.frame(res)
    sum_inc=sum(res$clin_inc0to80*1000)
    sum_inc_list[[j]] = sum_inc

  }
  check = unlist(sum_inc_list)
  check = as.data.frame(check)
  optimal = which.min(check$check)
  print(optimal)

  for (i in 1:500) {

    HR_profile = HR_data [i]
    HR_profile = as.data.frame(HR_profile)
    HR = HR_profile[,2]
    print (i)

    init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
    init_EIR <- 25

    time_period <- 365*5# run model for 5 years
    lower_time <- 365*3

    ivm_parms0 = ivm_fun(IVM_start_times = c(optimal+lower_time, optimal+30+lower_time, optimal+60+lower_time),  # time endectocide delivery occurs
                         time_period = time_period, # time period for the model to run over
                         hazard_profile = HR, # dummy hazard profile - must be vector
                         ivm_coverage = 0.8, # proportion of population receiving the endectocide
                         ivm_min_age = 5, # youngest age group receiving endectocide
                         ivm_max_age = 90) # oldest age group receiving endectocide

    wh <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                     num_int = 2,
                                     het_brackets = 5,
                                     age = init_age,
                                     init_EIR = 25,
                                     country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                     admin2 = "Fatick",
                                     ttt = ivm_parms0$ttt,
                                     eff_len = ivm_parms0$eff_len,
                                     haz = ivm_parms0$haz,
                                     ivm_cov_par = ivm_parms0$ivm_cov_par,
                                     ivm_min_age = ivm_parms0$ivm_min_age,
                                     ivm_max_age = ivm_parms0$ivm_max_age,
                                     IVRM_start = ivm_parms0$IVRM_start)

    res <- runfun(wh)
    res = as.data.frame(res)

    time = optimal+lower_time
    res$start_time= res$t-time

    res = res %>%
      select (t, start_time, clin_inc0to80, clin_inc0to5, clin_inc0to59, slide_prev0to80, slide_prev0to5, slide_prev0to59)

    data_frame_list [[i]] = res
  }
  try <- ldply(data_frame_list, data.frame)

  clin_inc0to80 = try %>% group_by (start_time) %>%  summarise(
    clin_05 = quantile(clin_inc0to80, probs=0.05, na.rm=TRUE),
    clin_25 = quantile(clin_inc0to80, probs=0.25, na.rm=TRUE),
    clin_med = quantile(clin_inc0to80, probs = 0.50, na.rm=TRUE),
    clin_75 = quantile(clin_inc0to80, probs = 0.75, na.rm=TRUE),
    clin_95 = quantile(clin_inc0to80, probs = 0.95, na.rm=TRUE))

  clin_inc0to5 = try %>% group_by (start_time) %>%  summarise(
    clin_05 = quantile(clin_inc0to5, probs=0.05, na.rm=TRUE),
    clin_25 = quantile(clin_inc0to5, probs=0.25, na.rm=TRUE),
    clin_med = quantile(clin_inc0to5, probs = 0.50, na.rm=TRUE),
    clin_75 = quantile(clin_inc0to5, probs = 0.75, na.rm=TRUE),
    clin_95 = quantile(clin_inc0to5, probs = 0.95, na.rm=TRUE))

  clin_inc0to59 = try %>% group_by (start_time) %>%  summarise(
    clin_05 = quantile(clin_inc0to59, probs=0.05, na.rm=TRUE),
    clin_25 = quantile(clin_inc0to59, probs=0.25, na.rm=TRUE),
    clin_med = quantile(clin_inc0to59, probs = 0.50, na.rm=TRUE),
    clin_75 = quantile(clin_inc0to59, probs = 0.75, na.rm=TRUE),
    clin_95 = quantile(clin_inc0to59, probs = 0.95, na.rm=TRUE))

  slide_prev0to80 = try %>% group_by (start_time) %>%  summarise(
    clin_05 = quantile(slide_prev0to80, probs=0.05, na.rm=TRUE),
    clin_25 = quantile(slide_prev0to80, probs=0.25, na.rm=TRUE),
    clin_med = quantile(slide_prev0to80, probs = 0.50, na.rm=TRUE),
    clin_75 = quantile(slide_prev0to80, probs = 0.75, na.rm=TRUE),
    clin_95 = quantile(slide_prev0to80, probs = 0.95, na.rm=TRUE))


  slide_prev0to5 = try %>% group_by (start_time) %>%  summarise(
    clin_05 = quantile(slide_prev0to5, probs=0.05, na.rm=TRUE),
    clin_25 = quantile(slide_prev0to5, probs=0.25, na.rm=TRUE),
    clin_med = quantile(slide_prev0to5, probs = 0.50, na.rm=TRUE),
    clin_75 = quantile(slide_prev0to5, probs = 0.75, na.rm=TRUE),
    clin_95 = quantile(slide_prev0to5, probs = 0.95, na.rm=TRUE))

  slide_prev0to59 = try %>% group_by (start_time) %>%  summarise(
    clin_05 = quantile(slide_prev0to59, probs=0.05, na.rm=TRUE),
    clin_25 = quantile(slide_prev0to59, probs=0.25, na.rm=TRUE),
    clin_med = quantile(slide_prev0to59, probs = 0.50, na.rm=TRUE),
    clin_75 = quantile(slide_prev0to59, probs = 0.75, na.rm=TRUE),
    clin_95 = quantile(slide_prev0to59, probs = 0.95, na.rm=TRUE))


  dataframes = list()

  dataframes[[1]] = clin_inc0to80
  dataframes[[2]] = clin_inc0to5
  dataframes[[3]] = clin_inc0to59
  dataframes[[4]] = slide_prev0to80
  dataframes[[5]] = slide_prev0to5
  dataframes[[6]] = slide_prev0to59


  return(dataframes)


}


IVM_1_400_vitro_DF = opt_date_Func (IVM_1_400_vitro)
IVM_3_300_vitro_DF = opt_date_Func (IVM_3_300_vitro)
IVM_3_600_vitro_DF = opt_date_Func (IVM_3_600_vitro)

IVM_1_400_vivo_DF = opt_date_Func (IVM_1_400_vivo)
IVM_3_300_vivo_DF = opt_date_Func (IVM_3_300_vivo)
IVM_3_600_vivo_DF = opt_date_Func (IVM_3_600_vivo)

NTBC_1_100_vitro_DF = opt_date_Func(NTBC_1_100_vitro)
NTBC_3_100_vitro_DF = opt_date_Func(NTBC_3_100_vitro)
NTBC_1_300_vitro_DF = opt_date_Func(NTBC_1_300_vitro)
NTBC_3_300_vitro_DF = opt_date_Func(NTBC_3_300_vitro)
NTBC_1_600_vitro_DF = opt_date_Func(NTBC_1_600_vitro)
NTBC_3_600_vitro_DF = opt_date_Func(NTBC_3_600_vitro)
NTBC_1_1mg_vitro_DF = opt_date_Func(NTBC_1_1mg_vitro)
NTBC_3_1mg_vitro_DF = opt_date_Func(NTBC_3_1mg_vitro)



DF_list<-lapply(list_all, FUN= opt_date_Func)


saveRDS(DF_list, file = "../01-Data/DF_list.RDS")


baseline_Func = function (data) {
  init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  init_EIR <- 25

  time_period <- 365*5# run model for 5 years
  HR = data [[1]]


  ivm_parms0 = ivm_fun(IVM_start_times = 365*8,  # time endectocide delivery occurs
                       time_period = time_period, # time period for the model to run over
                       hazard_profile = HR[,2], # dummy hazard profile - must be vector
                       ivm_coverage = 0.8, # proportion of population receiving the endectocide
                       ivm_min_age = 5, # youngest age group receiving endectocide
                       ivm_max_age = 90) # oldest age group receiving endectocide

  wh <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                   num_int = 2,
                                   het_brackets = 5,
                                   age = init_age,
                                   init_EIR = 25,
                                   country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                   admin2 = "Fatick",
                                   ttt = ivm_parms0$ttt,
                                   eff_len = ivm_parms0$eff_len,
                                   haz = ivm_parms0$haz,
                                   ivm_cov_par = ivm_parms0$ivm_cov_par,
                                   ivm_min_age = ivm_parms0$ivm_min_age,
                                   ivm_max_age = ivm_parms0$ivm_max_age,
                                   IVRM_start = ivm_parms0$IVRM_start)

  res <- runfun(wh)
  res = as.data.frame(res)

  res$start_time= res$t

  res = res %>%
    select (t, start_time, clin_inc0to80, clin_inc0to5, clin_inc0to59, slide_prev0to80, slide_prev0to5, slide_prev0to59)

  return (res)

}

baseline = baseline_Func (NTBC_1_100_vitro)


clin_list<-lapply(DF_list, FUN= clin_inc_Func)


clin_inc_Func = function (DF){

  for (i in 1:6) {

    data = (DF[[i]])

    baseline_DF  = baseline

    baseline_DF$start_time = data$start_time
    data = data %>%
      filter(start_time>0) %>%
      filter(start_time<366)
    data = data [,4]
    sum_inc = sum (data*1000)

    baseline_DF= baseline_DF %>%
      filter(start_time>0) %>%
      filter(start_time<366)

    baseline_DF = baseline_DF %>%
      select (-t, -start_time)

    baseline_DF = baseline_DF [,i]
    baseline_inc = sum (baseline_DF*1000)
    red = (baseline_inc - sum_inc)/ baseline_inc*100

    red_inc [[i]] = red
  }
  return (red_inc)

}



data = cbind (IVM_1_400_vt, IVM_3_300_vt, IVM_3_600_vt, IVM_1_400_vi, IVM_3_300_vi, IVM_3_600_vi,
              NTBC_1_100_vi, NTBC_3_100_vi, NTBC_1_300_vi, NTBC_3_300_vi, NTBC_1_600_vi, NTBC_3_600_vi, NTBC_1_1mg_vi, NTBC_3_1mg_vi)
rownames (data) = c("clin_inc0to80", "clin_inc0to5", "clin_inc0to59", "slide_prev0to80", "slide_prev0to5", "slide_prev0to59")

################
#Run transmission profiles 500 times


graph_Func = function (data, Title, outcome) {

  baseline_DF$start_time = data$start_time

  ggplot()+

    annotate("rect", xmin=0, xmax=1/30*3, ymin=0, ymax=2600 ,  alpha=0.8, fill="pink") +
    annotate("rect", xmin=1, xmax=1 + 1/30*3, ymin=0, ymax=2600 , alpha=0.8, fill="pink") +
    annotate("rect", xmin=2, xmax=2 + 1/30*3, ymin=0, ymax=2600 , alpha=0.8, fill="pink") +

    geom_line(data=baseline_DF, aes(x=start_time/30, y=outcome*1000*365, linetype="2"), size=0.7, col="grey") +
    geom_line(data=data, aes(x=start_time/30, y=clin_med*1000*365, linetype="1"), size=0.7) +
    geom_ribbon(data=data, aes(x=start_time/30, ymin = clin_05*1000*365, ymax = clin_95*1000*365, fill="25-75%"), alpha=0.4)+
    geom_ribbon(data=data, aes(x=start_time/30, ymin = clin_25*1000*365, ymax = clin_75*1000*365, fill="5-95%"), alpha=0.4)+
    xlab ("Months since intervention") +
    scale_x_continuous(expand = c(0,0), limits = c(-2, 20), breaks=seq(-2,20,2)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2600)) +
    ylab("All age clinical incidence")+
    scale_linetype_manual(values=c(1,2), labels = c("Median", "Baseline")) +
    #guides(linetype = guide_legend(override.aes = list(color = c("red", "red", "red")))) +
    scale_fill_manual (values=c('grey', 'grey'), labels=c("25-75%", "5-95%")) +

    ggtitle (Title) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size=9),
          axis.title.y = element_text(size=9),
          axis.title.x = element_text(size=9, vjust = -2),
          legend.title = element_blank(),legend.spacing.y = unit(0.1, 'mm'),
          axis.line = element_line(colour = "black", size = 0.1),
          legend.position = c(0.9,0.9),
          legend.background= element_rect(fill="white"),
          legend.text = element_text(size=9),
          plot.title = element_text(size=10, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour = "white"),
          legend.key = element_rect(fill = alpha("white", 0.0)),
          legend.text.align = 0,
          plot.margin = margin(5,5,5,5),
          legend.key.width = unit(0.3, "cm"))


}

IVM_1_400_vitro_DF
IVM_1_400_vitro_plot = graph_Func (IVM_1_400_vitro_DF[[1]], "IVM_1_400_vitro", baseline_DF$clin_inc0to5)

View(IVM_1_400_vitro_DF[[1]])

IVM_3_300_vitro_plot = graph_Func (IVM_3_300_vitro_DF[[1]], "IVM_3_300_vitro", baseline_DF$clin_inc0to5)


IVM_3_600_vitro_plot = graph_Func (IVM_3_600_vitro_DF, "IVM_3_600_vitro")


ggplot()+
  annotate("rect", xmin=0, xmax=1/30*3, ymin=0, ymax=2400 ,  alpha=0.8, fill="pink") +
  annotate("rect", xmin=1, xmax=1 + 1/30*3, ymin=0, ymax=2400 , alpha=0.8, fill="pink") +
  annotate("rect", xmin=2, xmax=2 + 1/30*3, ymin=0, ymax=2400 , alpha=0.8, fill="pink") +

  geom_line(data=res, aes(x=start_time/30, y=clin_inc0to80*1000*365, linetype="1"), color="red") +
  geom_ribbon(data=data_grouped, aes(x=start_time/30, ymin = clin_05*1000*365, ymax = clin_95*1000*365, fill="label_1"), alpha=0.5)+
  geom_ribbon(data=data_grouped, aes(x=start_time/30, ymin = clin_25*1000*365, ymax = clin_75*1000*365, fill="label_2"), alpha=0.5)+

  xlab ("Months since intervention") +
  scale_x_continuous(expand = c(0,0), limits = c(-2, 12), breaks=seq(-2,12,2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,2800)) +
  ylab("All age clinical incidence")+
  scale_linetype_manual(values=c(1), labels = c("Median")) +
  #guides(linetype = guide_legend(override.aes = list(color = c("red", "red", "red")))) +
  scale_fill_manual (values=c('blue', 'green'), labels=c("label_1", "label_2")) +

  ggtitle ("3X600 ug/kg IVM") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9, vjust = -2),
        legend.title = element_blank(),legend.spacing.y = unit(0.1, 'mm'),
        axis.line = element_line(colour = "black", size = 0.1),
        legend.position = c(0.83,0.9),
        legend.background= element_rect(fill="white"),
        legend.text = element_text(size=9),
        plot.title = element_text(size=10, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "white"),
        legend.key = element_rect(fill = alpha("white", 0.0)),
        legend.text.align = 0,
        plot.margin = margin(5,5,5,5),
        legend.key.width = unit(0.3, "cm"))



#generate dataframes per dosing regimen
seasonal_df=vector(mode="list",length=12)



# generate dataframes for baseline transmission
transmission_Func = function (start, cov, min, max, ivm_haz, eir, country, region) {

  ivm_parms = ivm_fun(IVM_start_times = start,
                      time_period = time_period,
                      hazard_profile = ivm_haz,
                      ivm_coverage=cov,
                      ivm_min_age=min,
                      ivm_max_age = max)

  wh <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                   num_int = 2,
                                   het_brackets = 5,
                                   age = init_age,
                                   init_EIR = eir,
                                   country = country,
                                   admin2 = region,
                                   ttt = ivm_parms$ttt,
                                   eff_len = ivm_parms$eff_len,
                                   haz = ivm_parms$haz,
                                   ivm_cov_par = ivm_parms$ivm_cov_par,
                                   ivm_min_age = ivm_parms$ivm_min_age,
                                   ivm_max_age = ivm_parms$ivm_max_age,
                                   IVRM_start = ivm_parms$IVRM_start)

  # Run all three scenarios
  res <- runfun(wh)
  res = as.data.frame(res)

  time = res$t

  res$date <- as.Date(time, origin = "2004-12-31")

  res$year <- as.numeric(format(res$date,'%Y'))

  res$month <- as.numeric(format(res$date,'%m'))



  res$time = res$t-3120

  res = res %>%
    group_by(year) %>%
    mutate (prev=mean(slide_prev0to80))

  return (res)

}

baseline_seasonal = transmission_Func (start = 20000, cov=0.7, min=5, max=95,
                                       ivm_haz= ivm_haz$NTBC_1_100, eir=25,
                                       country="Senegal", region="Fatick")

baseline_perennial = transmission_Func (start = 20000, cov=0.7, min=5, max=95,
                                        ivm_haz= ivm_haz$NTBC_1_100, eir=14,
                                        country="Democratic Republic of the Congo",
                                        region="Equateur")


new = baseline_perennial %>%
  select (year, prev)

new$prev = new$prev*100

new = distinct(new)

new = baseline_perennial %>%
  select (year, prev)

new$prev = new$prev*100

new = distinct(new)
#generate outcome parameters from simulations
incidence_sea=vector(mode="list",length=12)
prevalence_sea=vector(mode="list",length=12)

for  (i in 1:12) {

  data = seasonal_df [[i]]

  data = data %>%
    filter (t>3120) %>%
    filter (t<3485)

  sum_inc=sum(data$clin_inc0to80*1000)
  prev=mean(data$slide_prev0to80*1000)

  incidence_sea [i] = sum_inc
  prevalence_sea [i] = prev

}


incidence_per=vector(mode="list",length=12)
prevalence_per=vector(mode="list",length=12)

for  (i in 1:12) {

  data = perennial_df [[i]]

  data = data %>%
    filter (t>3120) %>%
    filter (t<3485)

  sum_inc=sum(data$clin_inc0to80*1000)
  prev=mean(data$slide_prev0to80*1000)

  incidence_per [i] = sum_inc
  prevalence_per [i] = prev

}



params_sea = cbind (incidence_sea, prevalence_sea, regimen)



params_sea$baseline_incidence = baseline_inc_sea
params_sea$baseline_prev = baseline_prev_sea
params_per$baseline_incidence = baseline_inc_per
params_per$baseline_prev = baseline_prev_per


params_per = params_per %>%
  mutate (inc_red = (baseline_inc_per-incidence_per)/baseline_inc_per*100) %>%
  mutate (prev_red = (baseline_prev_per - prevalence_per)/baseline_prev_per*100)

params_sea = params_sea %>%
  mutate (inc_red = (baseline_inc_sea-incidence_sea)/baseline_inc_sea*100) %>%
  mutate (prev_red = (baseline_prev_sea - prevalence_sea)/baseline_prev_sea*100)


params_sea = params_sea[ ,c(3,6,7)]
params_per = params_per[ ,c(3,6,7)]

View(data_grouped)
plot_Func = function (data_1, data_2, data_3, label_1, label_2, Title) {



  data_grouped$time = data_grouped$t-3122

  data_grouped = data_grouped %>%
    mutate(prev_05 = clin_05*1000*365) %>%
    mutate(prev_25 = clin_25*1000*365) %>%
    mutate(prev_50 = clin_med*1000*365) %>%
    mutate(prev_75 = clin_75*1000*365) %>%
    mutate(prev_95 = clin_95*1000*365)



  View(med_clin)
  plot =
    ggplot()+
    annotate("rect", xmin=0, xmax=1/30*3, ymin=0, ymax=2400 ,  alpha=0.8, fill="pink") +
    annotate("rect", xmin=1, xmax=1 + 1/30*3, ymin=0, ymax=2400 , alpha=0.8, fill="pink") +
    annotate("rect", xmin=2, xmax=2 + 1/30*3, ymin=0, ymax=2400 , alpha=0.8, fill="pink") +

    geom_line(data=res, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="1"), color="grey", size=0.8) +
    geom_line(data=res, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="2"), color="black", size=0.8) +
    geom_line(data=res, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="3"), color="grey", size=0.8) +

    xlab ("Months since intervention") +
    scale_x_continuous(expand = c(0,0), limits = c(-8, 12), breaks=seq(-2,20,2)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2800)) +
    ylab("All age clinical incidence")+
    scale_linetype_manual(values=c(1,1,1), labels = c("75% ci", "Median", "25% ci")) +
    guides(linetype = guide_legend(override.aes = list(color = c("grey", "black", "grey")))) +

    ggtitle ("3X600 ug/kg IVM") +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size=9),
          axis.title.y = element_text(size=9),
          axis.title.x = element_text(size=9, vjust = -2),
          legend.title = element_blank(),legend.spacing.y = unit(0.1, 'mm'),
          axis.line = element_line(colour = "black", size = 0.1),
          legend.position = c(0.83,0.9),
          legend.background= element_rect(fill="white"),
          legend.text = element_text(size=9),
          plot.title = element_text(size=10, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour = "white"),
          legend.key = element_rect(fill = alpha("white", 0.0)),
          legend.text.align = 0,
          plot.margin = margin(5,5,5,5),
          legend.key.width = unit(0.3, "cm"))
  return (plot)
}



plot_1 = plot_Func (seasonal_df[[9]], seasonal_df[[12]], baseline_seasonal, "IVM 3x300 (In-vitro)", "IVM 3x300 (In-vivo)", "Seasonal transmission")
plot_2 = plot_Func (perennial_df[[9]], perennial_df[[12]], baseline_perennial, "IVM 3x300 (In-vitro)", "IVM 3x300 (In-vivo)", "Perennial transmission")

plot_3 = plot_Func (seasonal_df[[4]], seasonal_df[[9]], baseline_seasonal, "NTBC 3x300", "IVM 3x300", "Seasonal transmission")
plot_4 = plot_Func (perennial_df[[4]], perennial_df[[9]], baseline_perennial, "NTBC 3x300", "IVM 3x300", "Perennial transmission")

plot_5 = plot_Func (seasonal_df[[3]], seasonal_df[[4]], baseline_seasonal, "NTBC 1x300", "NTBC 3x300", "Seasonal transmission")
plot_6 = plot_Func (perennial_df[[3]], perennial_df[[4]], baseline_perennial, "NTBC 1x300", "NTBC 3x300", "Perennial transmission")

plot_7 = plot_Func (seasonal_df[[1]], seasonal_df[[2]], baseline_seasonal, "NTBC 1x100", "NTBC 3x100", "Seasonal transmission")
plot_8 = plot_Func (perennial_df[[1]], perennial_df[[2]], baseline_perennial, "NTBC 1x100", "NTBC 3x100", "Perennial transmission")


#Make time series

res0 = as.data.frame(res0)

start = as.Date("2010-01-01")
res0$date = as.Date(res0$t, origin=start)




