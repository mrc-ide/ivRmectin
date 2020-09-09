

# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)


# Create a vector of age categories for the model
init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 3

#low - 2, moderate - 15, high - 120


# Provide the length of time (in days) that you want to run the model for
time_period <- 3650 # run model for 10 years

# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

# Load ivermectin hazardd
ivm_haz_1 <- read.table("data/ivermectin_hazards.txt", header=TRUE)
ivm_haz_2 <- read.csv("data/sim_ave_HR.csv", header=TRUE)
ivm_haz_2 = ivm_haz_2[ ,2:6]

ivm_haz_1 = ivm_haz_1[1:28,]
ivm_haz_2 = ivm_haz_2[1:28, ]
ivm_haz=cbind(ivm_haz_1, ivm_haz_2)

colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS","IVM_300_3_AT", "IVM_600_3_AT", "IVM_400_1_AT", "NTBC_1000_1_AT", "NTBC_1000_3_AT")


# Running the ivm_fun function to generate the extra endectocide specific parameters that
# you have to pass to the model
ivm_parms = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                     time_period = time_period, # time period for the model to run over
                     hazard_profile = ivm_haz$IVM_300_3_AT[1:23], # dummy hazard profile - must be vector
                     ivm_coverage = 0.8, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide

# Creates the odin model with all the required parameters - it is then ready to run
# although this isn't the part where the model is actually run - that's below.
wh <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                num_int = 4, # number of vector control (IRS and ITN) population groups
                                het_brackets = 5, # number of heterogeneous biting categories
                                age = init_age, # the different age classes to be ran within the model
                                init_EIR = init_EIR, # the Entomological Innoculation Rate
                                country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                ttt = ivm_parms$ttt, # model specific parameter to control timing of endectocide delivery
                                eff_len = ivm_parms$eff_len, # number of days after receiving endectocide that HR is higher
                                haz = ivm_parms$haz, # hazard ratio for each off the eff_len number of days
                                ivm_cov_par = ivm_parms$ivm_cov_par, # proportion of popuulation receiving the endectocide
                                ivm_min_age = ivm_parms$ivm_min_age, # youngest age group receiving endectocide
                                ivm_max_age = ivm_parms$ivm_max_age, # oldest age group receiving endectocide
                                IVRM_start = ivm_parms$IVRM_start) # model specific parameter to control timing of endectocide delivery

# This function runs the odin model with all the parameters that you previously generated with
# the ivRmectimn::create_r_model call
runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

res0 <- runfun(wh)


#loop to obtain optimal timing of seasonal MDA
#warning - takes days to run
##Function - Optimisation

vector=vector(mode="numeric",length=3)
list=vector(mode="list",length=3650)


for (i in 1:3650) {
  vector[1] <- i
  vector[2] <- vector[1] + 30
  vector[3] <- vector[2] + 30
  list[[i]] <- vector
}

sum_inc_list=vector(mode="list",length=3650)


for (i in 3000:3650) {
    # Create a vector of age categories for the model
  init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  # Provide a value of the annual EIR for this model run
  init_EIR <- 14
  #low - 2, moderate - 15, high - 120
  # Provide the length of time (in days) that you want to run the model for
  time_period <- 365 *10 # run model for 10 years
  vector<-list[i]

 administration = unlist(vector)
   ivm_parms0 = ivm_fun(IVM_start_times = administration,  # time endectocide delivery occurs
                       time_period = time_period, # time period for the model to run over
                       hazard_profile = ivm_haz$IVM_300_3_AT, # dummy hazard profile - must be vector
                       ivm_coverage = 0.8, # proportion of population receiving the endectocide
                       ivm_min_age = 5, # youngest age group receiving endectocide
                       ivm_max_age = 90) # oldest age group receiving endectocide

  wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 14,
                                    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    admin2 = "Fatick",
                                    ttt = ivm_parms0$ttt,
                                    eff_len = ivm_parms0$eff_len,
                                    haz = ivm_parms0$haz,
                                    ivm_cov_par = ivm_parms0$ivm_cov_par,
                                    ivm_min_age = ivm_parms0$ivm_min_age,
                                    ivm_max_age = ivm_parms0$ivm_max_age,
                                    IVRM_start = ivm_parms0$IVRM_start)

  res0 <- runfun(wh)

  sum_inc=sum(res$clin_inc0to80*1000*365)

  sum_inc_list[i] <- sum_inc

}

#Obtain optimal start dates combination
sum_inc_list = administration

administration = unlist(sum_inc_list)
View(administration)
min(administration)


#Make time series

res0 = as.data.frame(res0)

start = as.Date("2010-01-01")
res0$date = as.Date(res0$t, origin=start)

#####Functions for generating graphs and output parameters#################
#Function for generating clinical incidence and slide prevalence reductions for three scenarios
Params_Func = function (eir, start1, start2, Country, Region,
                        HR_sim0, cov0, min0, max0,
                        HR_sim1, cov1, min1, max1,
                        HR_sim2, cov2, min2, max2) {

  # no ivermectin being used
  ivm_parms0 = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                       time_period = time_period, # time period for the model to run over
                       hazard_profile = HR_sim0, # dummy hazard profile - must be vector
                       ivm_coverage = cov0, # proportion of population receiving the endectocide
                       ivm_min_age = min0, # youngest age group receiving endectocide
                       ivm_max_age = max0) # oldest age group receiving endectocide


  wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    admin2 = Region,
                                    ttt = ivm_parms0$ttt,
                                    eff_len = ivm_parms0$eff_len,
                                    haz = ivm_parms0$haz,
                                    ivm_cov_par = ivm_parms0$ivm_cov_par,
                                    ivm_min_age = ivm_parms0$ivm_min_age,
                                    ivm_max_age = ivm_parms0$ivm_max_age,
                                    IVRM_start = ivm_parms0$IVRM_start)

  res <- runfun(wh)


  # First endectocide scenario
  ivm_parms1 = ivm_fun(IVM_start_times = start1,
                       time_period = time_period,
                       hazard_profile = HR_sim1,
                       ivm_coverage=cov1,
                       ivm_min_age=min1,
                       ivm_max_age = max1)


  wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start)

  # Second endectocide scenario
  ivm_parms2 = ivm_fun(IVM_start_times = start2,
                       time_period = time_period,
                       hazard_profile = HR_sim2,
                       ivm_coverage=cov2,
                       ivm_min_age=min2,
                       ivm_max_age = max2)

  wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms2$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start)

  # Run all three scenarios
  res0 <- runfun(wh0)
  res1 <- runfun(wh1)
  res2 <- runfun(wh2)

    res0 = as.data.frame(res0)
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)

  res0 = res0 %>%
     filter (t>3120) %>%
    filter (t<3485)

  res1 = res1 %>%
    filter (t>3120) %>%
    filter (t<3485)

  res2 = res2 %>%
    filter (t>3120) %>%
    filter (t<3485)

  # Plotting the results
  sum_inc0=sum(res0$clin_inc0to80*1000)
  sum_inc1=sum(res1$clin_inc0to80*1000)
  sum_inc2=sum(res2$clin_inc0to80*1000)

  inc_red1= (sum_inc0 - sum_inc1)/sum_inc0*100
  inc_red2= (sum_inc0 - sum_inc2)/sum_inc0*100

  ave_prev0=mean(res0$slide_prev0to80*100)
  ave_prev1=mean(res1$slide_prev0to80*100)
  ave_prev2=mean(res2$slide_prev0to80*100)

  prev_red1= (ave_prev0 - ave_prev1)/ave_prev0*100
  prev_red2= (ave_prev0 - ave_prev2)/ave_prev0*100

  outcome_params <- c(sum_inc0, sum_inc1, sum_inc2, inc_red1, inc_red2, ave_prev0, ave_prev1, ave_prev2,  prev_red1, prev_red2)

  return (outcome_params)

}



#Function which returns total incidence and prevalence statistics for comparison
Inc_plot_Func_2 = function (eir, start1, start2, Country, Region,
                            HR_sim0, cov0, min0, max0,
                            HR_sim1, cov1, min1, max1,
                            HR_sim2, cov2, min2, max2,
                            ylim1, ylim2, labels, Title, xlablabel, ylablabel) {

  # no ivermectin being used
  ivm_parms0 = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                       time_period = time_period, # time period for the model to run over
                       hazard_profile = HR_sim0, # dummy hazard profile - must be vector
                       ivm_coverage = cov0, # proportion of population receiving the endectocide
                       ivm_min_age = min0, # youngest age group receiving endectocide
                       ivm_max_age = max0) # oldest age group receiving endectocide


  wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    admin2 = Region,
                                    ttt = ivm_parms0$ttt,
                                    eff_len = ivm_parms0$eff_len,
                                    haz = ivm_parms0$haz,
                                    ivm_cov_par = ivm_parms0$ivm_cov_par,
                                    ivm_min_age = ivm_parms0$ivm_min_age,
                                    ivm_max_age = ivm_parms0$ivm_max_age,
                                    IVRM_start = ivm_parms0$IVRM_start)

  # First endectocide scenario
  ivm_parms1 = ivm_fun(IVM_start_times = start1,
                       time_period = time_period,
                       hazard_profile = HR_sim1,
                       ivm_coverage=cov1,
                       ivm_min_age=min1,
                       ivm_max_age = max1)


  wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start)

  # Second endectocide scenario
  ivm_parms2 = ivm_fun(IVM_start_times = start2,
                       time_period = time_period,
                       hazard_profile = HR_sim2,
                       ivm_coverage=cov2,
                       ivm_min_age=min2,
                       ivm_max_age = max2)

  wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms2$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start)

  # Run all three scenarios
  res0 <- runfun(wh0)
  res1 <- runfun(wh1)
  res2 <- runfun(wh2)


  res0$time = res0$t-3120
  res1$time = res1$t-3120
  res2$time = res2$t-3120


  res0 = as.data.frame(res0)
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)



  plot =  ggplot()+

    geom_line(data=res1, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="1"), color="#689ed4", size=0.2) +
    geom_line(data=res2, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="2"), color="#c07142", size=0.2) +

    xlab (xlablabel) +
    scale_x_continuous(expand = c(0,0), limits = c(-2,13), breaks=seq(-2,13,2)) +
    scale_y_continuous(expand = c(0,0), limits = ylim1) +
    ylab(ylablabel)+
    scale_linetype_manual(values=c(1,1), labels = labels) +
    guides(linetype = guide_legend(override.aes = list(color = c("#689ed4", "#c07142")))) +

    annotate("rect", xmin=0, xmax=1/30*3, ymin=0, ymax=1300 ,  alpha=0.2, col="pink", fill="pink") +
    annotate("rect", xmin=1, xmax=1 + 1/30*3, ymin=0, ymax=1300 , alpha=0.2, col="pink", fill="pink") +
    annotate("rect", xmin=2, xmax=2 + 1/30*3, ymin=0, ymax=1300 , alpha=0.2, col="pink", fill="pink") +

    ggtitle (Title) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(size = 4),
          axis.text.y = element_text(size=4),
          axis.title.y = element_text(size=4),
          axis.title.x = element_text(size=4, vjust = -2),
          legend.title = element_blank(),legend.spacing.y = unit(0.1, 'mm'),
          axis.line = element_line(colour = "black", size = 0.1),
          legend.position = c(0.7,0.9),
          legend.background= element_rect(fill="white"),
          legend.text = element_text(size=3),
          plot.title = element_text(size=6, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour = "white"),
          legend.key = element_rect(fill = alpha("white", 0.0)),
          legend.text.align = 0,
          plot.margin = margin(5,5,5,5),
          legend.key.width = unit(0.5, "cm"))


  return (plot)

}





#Function which returns total incidence and prevalence statistics for comparison
Inc_plot_Func_3 = function (eir, start1, start2, Country, Region,
                        HR_sim0, cov0, min0, max0,
                        HR_sim1, cov1, min1, max1,
                        HR_sim2, cov2, min2, max2,
                        ylim1, ylim2, labels, Title, xlablabel, ylablabel) {

  # no ivermectin being used
  ivm_parms0 = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                       time_period = time_period, # time period for the model to run over
                       hazard_profile = HR_sim0, # dummy hazard profile - must be vector
                       ivm_coverage = cov0, # proportion of population receiving the endectocide
                       ivm_min_age = min0, # youngest age group receiving endectocide
                       ivm_max_age = max0) # oldest age group receiving endectocide


  wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    admin2 = Region,
                                    ttt = ivm_parms0$ttt,
                                    eff_len = ivm_parms0$eff_len,
                                    haz = ivm_parms0$haz,
                                    ivm_cov_par = ivm_parms0$ivm_cov_par,
                                    ivm_min_age = ivm_parms0$ivm_min_age,
                                    ivm_max_age = ivm_parms0$ivm_max_age,
                                    IVRM_start = ivm_parms0$IVRM_start)


  # First endectocide scenario
  ivm_parms1 = ivm_fun(IVM_start_times = start1,
                       time_period = time_period,
                       hazard_profile = HR_sim1,
                       ivm_coverage=cov1,
                       ivm_min_age=min1,
                       ivm_max_age = max1)


  wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start)

  # Second endectocide scenario
  ivm_parms2 = ivm_fun(IVM_start_times = start2,
                       time_period = time_period,
                       hazard_profile = HR_sim2,
                       ivm_coverage=cov2,
                       ivm_min_age=min2,
                       ivm_max_age = max2)

  wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms2$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start)

  # Run all three scenarios
  res0 <- runfun(wh0)
  res1 <- runfun(wh1)
  res2 <- runfun(wh2)


  res0 = as.data.frame(res0)
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)

  res0$time = res0$t-3120
  res1$time = res1$t-3120
  res2$time = res2$t-3120



  plot =  ggplot()+

    geom_line(data=res0, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="3"), color="darkgrey", size=0.2) +

    geom_line(data=res1, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="1"), color="#689ed4", size=0.2) +
    geom_line(data=res2, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="2"), color="#c07142", size=0.2) +

    xlab (xlablabel) +
    scale_x_continuous(expand = c(0,0), limits = c(-2,13), breaks=seq(-2,13,2)) +
    scale_y_continuous(expand = c(0,0), limits = ylim1) +
    ylab(ylablabel)+
    scale_linetype_manual(values=c(1,1,2), labels = labels) +
    guides(linetype = guide_legend(override.aes = list(color = c("#689ed4", "#c07142", "darkgrey")))) +

    annotate("rect", xmin=0, xmax=1/30*3, ymin=0, ymax=2300 ,  alpha=0.2, col="pink", fill="pink") +
    annotate("rect", xmin=1, xmax=1 + 1/30*3, ymin=0, ymax=2300 , alpha=0.2, col="pink", fill="pink") +
    annotate("rect", xmin=2, xmax=2 + 1/30*3, ymin=0, ymax=2300 , alpha=0.2, col="pink", fill="pink") +

    ggtitle (Title) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(size = 4),
          axis.text.y = element_text(size=4),
          axis.title.y = element_text(size=4),
          axis.title.x = element_text(size=4, vjust = -2),
          legend.title = element_blank(),legend.spacing.y = unit(0.1, 'mm'),
          axis.line = element_line(colour = "black", size = 0.1),
          legend.position = c(0.7,0.9),
          legend.background= element_rect(fill="white"),
          legend.text = element_text(size=3),
          plot.title = element_text(size=6, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour = "white"),
          legend.key = element_rect(fill = alpha("white", 0.0)),
          legend.text.align = 0,
          plot.margin = margin(5,5,5,5),
          legend.key.width = unit(0.5, "cm"))


  return (plot)

}



Prev_plot_Func_2 = function (eir, start1, start2, Country, Region,
                             HR_sim0, cov0, min0, max0,
                             HR_sim1, cov1, min1, max1,
                             HR_sim2, cov2, min2, max2,
                             ylim1, ylim2, labels, Title) {

  # no ivermectin being used
  ivm_parms0 = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                       time_period = time_period, # time period for the model to run over
                       hazard_profile = HR_sim0, # dummy hazard profile - must be vector
                       ivm_coverage = cov0, # proportion of population receiving the endectocide
                       ivm_min_age = min0, # youngest age group receiving endectocide
                       ivm_max_age = max0) # oldest age group receiving endectocide


  wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    admin2 = Region,
                                    ttt = ivm_parms0$ttt,
                                    eff_len = ivm_parms0$eff_len,
                                    haz = ivm_parms0$haz,
                                    ivm_cov_par = ivm_parms0$ivm_cov_par,
                                    ivm_min_age = ivm_parms0$ivm_min_age,
                                    ivm_max_age = ivm_parms0$ivm_max_age,
                                    IVRM_start = ivm_parms0$IVRM_start)

  res0 <- runfun(wh)


  # First endectocide scenario
  ivm_parms1 = ivm_fun(IVM_start_times = start1,
                       time_period = time_period,
                       hazard_profile = HR_sim1,
                       ivm_coverage=cov1,
                       ivm_min_age=min1,
                       ivm_max_age = max1)


  wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start)

  # Second endectocide scenario
  ivm_parms2 = ivm_fun(IVM_start_times = start2,
                       time_period = time_period,
                       hazard_profile = HR_sim2,
                       ivm_coverage=cov2,
                       ivm_min_age=min2,
                       ivm_max_age = max2)

  wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms2$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start)

  # Run all three scenarios
  res0 <- runfun(wh0)
  res1 <- runfun(wh1)
  res2 <- runfun(wh2)


  res0$time = res0$t-3120
  res1$time = res1$t-3120
  res2$time = res2$t-3120

  res0 = as.data.frame(res0)
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)


  plot =
    ggplot()+

    geom_line(data=res1, aes(x=time/30, y=slide_prev0to80*100, linetype="1"), color="#689ed4", size=0.2) +
    geom_line(data=res2, aes(x=time/30, y=slide_prev0to80*100, linetype="2"), color="#c07142", size=0.2) +

    xlab ("Time since start of intervention (Months)") +
    scale_x_continuous(expand = c(0,0), limits = c(-2,13), breaks=seq(-2,13,2)) +
    scale_y_continuous(expand = c(0,0), limits = ylim2) +
    ylab(ylablabel)+
    scale_linetype_manual(values=c(1, 1), labels = labels) +
    guides(linetype = guide_legend(override.aes = list(color = c("#689ed4", "#c07142")))) +

    annotate("rect", xmin=0 , xmax=1/30*3, ymin=0, ymax=50 ,   col="pink", fill="pink", alpha=0.2) +
    annotate("rect", xmin=1, xmax=1 + 1/30*3, ymin=0, ymax=50 , col="pink", fill="pink", alpha=0.2) +
    annotate("rect", xmin=2, xmax=2 + 1/30*3, ymin=0, ymax=50 , col="pink", fill="pink", alpha=0.2) +
    ggtitle (Title) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(size = 4),
          axis.text.y = element_text(size=4),
          axis.title.y = element_text(size=4),
          axis.title.x = element_text(size=4, vjust = -2),
          legend.title = element_blank(),legend.spacing.y = unit(0.1, 'mm'),
          axis.line = element_line(colour = "black", size = 0.1),
          legend.position = c(0.7,0.9),
          legend.background= element_rect(fill="white"),
          legend.text = element_text(size=3),
          plot.title = element_text(size=3, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour = "white"),
          legend.key = element_rect(fill = alpha("white", 0.0)),
          legend.text.align = 0,
          plot.margin = margin(5,5,5,5),
          legend.key.width = unit(0.5, "cm"))



  return (plot)

}


Prev_plot_Func_3 = function (eir, start1, start2, Country, Region,
                             HR_sim0, cov0, min0, max0,
                             HR_sim1, cov1, min1, max1,
                             HR_sim2, cov2, min2, max2,
                             ylim1, ylim2, labels, Title, ylablabel) {

  # no ivermectin being used
  ivm_parms0 = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                       time_period = time_period, # time period for the model to run over
                       hazard_profile = HR_sim0, # dummy hazard profile - must be vector
                       ivm_coverage = cov0, # proportion of population receiving the endectocide
                       ivm_min_age = min0, # youngest age group receiving endectocide
                       ivm_max_age = max0) # oldest age group receiving endectocide


  wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    admin2 = Region,
                                    ttt = ivm_parms0$ttt,
                                    eff_len = ivm_parms0$eff_len,
                                    haz = ivm_parms0$haz,
                                    ivm_cov_par = ivm_parms0$ivm_cov_par,
                                    ivm_min_age = ivm_parms0$ivm_min_age,
                                    ivm_max_age = ivm_parms0$ivm_max_age,
                                    IVRM_start = ivm_parms0$IVRM_start)

  res0 <- runfun(wh)


  # First endectocide scenario
  ivm_parms1 = ivm_fun(IVM_start_times = start1,
                       time_period = time_period,
                       hazard_profile = HR_sim1,
                       ivm_coverage=cov1,
                       ivm_min_age=min1,
                       ivm_max_age = max1)


  wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start)

  # Second endectocide scenario
  ivm_parms2 = ivm_fun(IVM_start_times = start2,
                       time_period = time_period,
                       hazard_profile = HR_sim2,
                       ivm_coverage=cov2,
                       ivm_min_age=min2,
                       ivm_max_age = max2)

  wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 2,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = eir,
                                    country = Country,
                                    admin2 = Region,
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms2$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start)

  # Run all three scenarios
  res0 <- runfun(wh0)
  res1 <- runfun(wh1)
  res2 <- runfun(wh2)


  res0$time = res0$t-3120
  res1$time = res1$t-3120
  res2$time = res2$t-3120

  res0 = as.data.frame(res0)
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)


  plot =
    ggplot()+

    geom_line(data=res0, aes(x=time/30, y=slide_prev0to80*100, linetype="3"), color="darkgrey", size=0.2) +

    geom_line(data=res1, aes(x=time/30, y=slide_prev0to80*100, linetype="1"), color="#689ed4", size=0.2) +
    geom_line(data=res2, aes(x=time/30, y=slide_prev0to80*100, linetype="2"), color="#c07142", size=0.2) +

    xlab ("Time since start of intervention (Months)") +
    scale_x_continuous(expand = c(0,0), limits = c(-2,13), breaks=seq(-2,13,2)) +
    scale_y_continuous(expand = c(0,0), limits = ylim2) +
    ylab(ylablabel)+
    scale_linetype_manual(values=c(1, 1, 2), labels = labels) +
    guides(linetype = guide_legend(override.aes = list(color = c("#689ed4", "#c07142", "darkgrey")))) +

    annotate("rect", xmin=0 , xmax=1/30*3, ymin=0, ymax=50 ,   col="pink", fill="pink", alpha=0.2) +
    annotate("rect", xmin=1, xmax=1 + 1/30*3, ymin=0, ymax=50 , col="pink", fill="pink", alpha=0.2) +
    annotate("rect", xmin=2, xmax=2 + 1/30*3, ymin=0, ymax=50 , col="pink", fill="pink", alpha=0.2) +
    ggtitle (Title) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(size = 4),
          axis.text.y = element_text(size=4),
          axis.title.y = element_text(size=4),
          axis.title.x = element_text(size=4, vjust = -2),
          legend.title = element_blank(),legend.spacing.y = unit(0.1, 'mm'),
          axis.line = element_line(colour = "black", size = 0.1),
          legend.position = c(0.7,0.9),
          legend.background= element_rect(fill="white"),
          legend.text = element_text(size=3),
          plot.title = element_text(size=3, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour = "white"),
          legend.key = element_rect(fill = alpha("white", 0.0)),
          legend.text.align = 0,
          plot.margin = margin(5,5,5,5),
          legend.key.width = unit(0.5, "cm"))


  return (plot)

}

#graphical output

View(ivm_haz)
seasonal_1 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov2=0.7, min2=5, max2=90)


perennial_1 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov2=0.7, min2=5, max2=90)



plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov2=0.7, min2=5, max2=90,
                         ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 1x400ug/kg (Slater et al)", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                         Title= "Highly seasonal transmission", xlablabel="", ylablabel="Clinical incidence per 1,000 population")

plot_2 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                          HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 1x400ug/kg (Slater et al)", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                          Title= "Perennial  transmission", xlablabel="", ylablabel="")


 plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov2=0.7, min2=5, max2=90,
                         ylim1=c(0, 1000), ylim2=c(0,50), labels=c("IVM 1x400ug/kg", "IVM 3x300ug/kg"),
                         Title= "Highly seasonal transmission", xlablabel= "Time since start of intervention (Months)", ylablabel= "Annual slie prevalence (%)")

 plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 1000), ylim2=c(0,50), labels=c("IVM 1x400ug/kg", "IVM 3x300ug/kg"),
                           Title= "Perennial transmission", xlablabel= "Time since start of intervention (Months)", ylablabel="")



 plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                          HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 1x400ug/kg (Slater et al)", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                          Title= "", ylablabel="Annual slide prevalence")

 plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                          HR_sim0=ivm_haz$IVM_300_3_AT_3_HS[1:23],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_400_1_AT_1_HS[1:23], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:23], cov2=0.7, min2=0, max2=90,
                          ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 1x400ug/kg (Slater et al)", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                          Title= "", ylablabel= "" )

 grid_1 = grid.arrange(plot_1, plot_2, plot_5, plot_6, nrow = 2, ncol=2)

 ggsave(file="grid_1.png", grid_1, height=4, width=4, dpi=400, limitsize=FALSE)

 grid_1a = grid.arrange(plot_3, plot_4, nrow = 1, ncol=2)
 ggsave(file="grid_1a.png", grid_1a, height=2, width=3, dpi=400, limitsize=FALSE)

dev.off()
####################################

seasonal_2 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_600[1:28], cov2=0.7, min2=5, max2=90)


perennial_2 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_600[1:28], cov2=0.7, min2=5, max2=90)





plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                          Title= "Highly seasonal transmission", xlablabel="", ylablabel="Clinical incidence per 1,000 population")

plot_2 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                          Title= "Perennial  transmission", xlablabel="", ylablabel="")


plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 1300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)"),
                          Title= "Highly seasonal transmission", xlablabel="Time since start of intervention (Months)", ylablabel="Clinical incidence per 1,000 population")

plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 1300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)"),
                          Title= "Perennial transmission", xlablabel="Time since start of intervention (Months)", ylablabel="")

plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                           Title= "", ylablabel="Annual slide prevalence")


plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                           Title= "", ylablabel= "" )

grid_2 = grid.arrange(plot_1, plot_2, plot_5, plot_6, nrow = 2, ncol=2)

ggsave(file="grid_2.png", grid_2, height=4, width=4, dpi=400, limitsize=FALSE)


grid_2a = grid.arrange(plot_3, plot_4, nrow = 1, ncol=2)

ggsave(file="grid_2a.png", grid_2a, height=2, width=3.5, dpi=400, limitsize=FALSE)

grid_2b = grid.arrange(plot_5, plot_6, nrow = 1, ncol=2)

ggsave(file="grid_2b.png", grid_2b, height=2, width=3.5, dpi=400, limitsize=FALSE)

###################

seasonal_3 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                        HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                        HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                        HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90)


perennial_3 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                        HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                        HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                        HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90)



plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                          Title= "Highly seasonal transmission", xlablabel="", ylablabel="Clinical incidence per 1,000 population")

plot_2 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                          Title= "Perennial  transmission", xlablabel="", ylablabel="")


plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 1300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)"),
                          Title= "Highly seasonal transmission")

plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                          HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                          HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                          HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                          ylim1=c(0, 1300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)"),
                          Title= "Perennial transmission")

plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                           Title= "", ylablabel="Annual slide prevalence")


plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                           Title= "", ylablabel="")

grid_3 = grid.arrange(plot_1, plot_2, plot_5, plot_6, nrow = 2, ncol=2)

ggsave(file="grid_3.png", grid_3, height=4, width=4, dpi=400, limitsize=FALSE)


#########
###################


 ##########

seasonal_4 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90)


perennial_4 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90)




 plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                           Title= "Highly seasonal transmission", xlablabel="", ylablabel="Clinical incidence per 1,000 population")

 plot_2 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("IVM 300ug/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                           Title= "Perennial  transmission", xlablabel="", ylablabel="")


 plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 1300), ylim2=c(0,50), labels=c("IVM 3x300ug/kg","NTBC 3x1mg/kg"),
                           Title= "Highly seasonal transmission")

 plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "NTBC 3x1mg/kg"),
                           Title= "Perennial  transmission")


 plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                            HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                            Title= "", ylablabel="Annual slide prevalence")


 plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                            HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$IVM_300_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=5, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                            Title= "", ylablabel="")

 grid_4 = grid.arrange(plot_1, plot_2, plot_5, plot_6, nrow = 2, ncol=2)

 ggsave(file="grid_4.png", grid_4, height=4, width=4, dpi=400, limitsize=FALSE)
 #################

 seasonal_5 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90)


 perennial_5 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                            HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90)



 plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                           Title= "Highly seasonal transmission", xlablabel="", ylablabel="Clinical incidence per 1,000 population")


 plot_2 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                           Title= "Perennial  transmission", xlablabel="", ylablabel="")


 plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 1300), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "IVM 3x300ug/kg (Slater et al)",),
                           Title= "Highly seasonal transmission")

 plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "IVM 3x300ug/kg (Slater et al)",),
                           Title= "Perennial  transmission")


 plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                            HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                            Title= "", ylablabel="Annual slide prevalence")


 plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                            HR_sim0=ivm_haz$IVM_300_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$IVM_300_3_AT_3_HS[1:28], cov2=0.7, min2=5, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "IVM 3x300ug/kg (Slater et al)", "Baseline malaria"),
                            Title= "", ylablabel="")

 grid_5 = grid.arrange(plot_1, plot_2, plot_5, plot_6, nrow = 2, ncol=2)

 ggsave(file="grid_5.png", grid_5, height=4, width=4, dpi=400, limitsize=FALSE)

 ##########

 seasonal_6 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90)


 perennial_6 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90)

 df = rbind(seasonal_1, perennial_1, seasonal_2, perennial_2, seasonal_3, perennial_3, seasonal_4, perennial_4, seasonal_5, perennial_5, seasonal_6, perennial_6)
 df = as.data.frame(df)
 colnames (df) = c("sum_inc0", "sum_inc1", "sum_inc2", "inc_red1", "inc_red2", "ave_prev0", "ave_prev1", "ave_prev2",  "prev_red1", "prev_red2")
View(df)


 plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("NTBC (5-90yrs)", "NTBC (0-90yrs)", "Baseline malaria"),
                           Title= "Highly seasonal transmission", xlablabel="", ylablabel="Clinical incidence per 1,000 population")

 plot_2 = Inc_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("NTBC (5-90yrs)", "NTBC (0-90yrs)", "Baseline malaria"),
                           Title= "Perennial  transmission",  xlablabel="", ylablabel="")


 plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("NTBC (5-90yrs)", "NTBC (0-90yrs)"),
                           Title= "Highly seasonal transmission")

 plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("NTBC (5-90yrs)", "NTBC (0-90yrs)"),
                           Title= "Perennial  transmission")


 plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                            HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("NTBC (5-90yrs)", "NTBC (0-90yrs)", "Baseline malaria"),
                            Title= "", ylablabel="Annual slide prevalence")


 plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                            HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_3_AT[1:28], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("NTBC (5-90yrs)", "NTBC (0-90yrs)", "Baseline malaria"),
                            Title= "", ylablabel="")

##################


 ##########

 seasonal_7 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                           HR_sim2=ivm_haz$IVM_400_1_AT[1:28], cov2=0.7, min2=5, max2=90)


 perennial_7 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                            HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                            HR_sim2=ivm_haz$IVM_400_1_AT[1:28], cov2=0.7, min2=5, max2=90)

 df = rbind(seasonal_7, perennial_7)
 , seasonal_2, perennial_2, seasonal_3, perennial_3, seasonal_4, perennial_4, seasonal_5, perennial_5, seasonal_6, perennial_6)
 df = as.data.frame(df)
 colnames (df) = c("sum_inc0", "sum_inc1", "sum_inc2", "inc_red1", "inc_red2", "ave_prev0", "ave_prev1", "ave_prev2",  "prev_red1", "prev_red2")
 View(df)


 plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                           Title= "Highly seasonal transmission", xlablabel="", ylablabel="Clinical incidence per 1,000 population")

 plot_2 = Inc_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 2300), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                           Title= "Perennial  transmission",  xlablabel="", ylablabel="")


 plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "NTBC 3x1mg/kg",),
                           Title= "Highly seasonal transmission")

 plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "NTBC 3x1mg/kg",),
                           Title= "Perennial  transmission")


 plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                            HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                            Title= "", ylablabel="Annual slide prevalence")


 plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                            HR_sim0=ivm_haz$NTBC_1000_3_AT[1:28],cov0=0.7, min0=0, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_1_AT[1:28], cov1=0.7, min1=0, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:28], cov2=0.7, min2=0, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("NTBC 1x1mg/kg", "NTBC 3x1mg/kg", "Baseline malaria"),
                            Title= "", ylablabel="")


 grid_7 = grid.arrange(plot_1, plot_2, plot_5, plot_6, nrow = 2, ncol=2)

 ggsave(file="grid_7.png", grid_7, height=4, width=4, dpi=400, limitsize=FALSE)




 #Hypothetical endectocides:
 hr_0= rep(0,28)

 hr_1= rep(1,28)
 hr_2= rep(2,28)
 hr_3= rep(3,28)
 hr_4= rep(4,28)
 hr_5= rep(5,28)
 hr_6= rep(6,28)
 hr_7= rep(7,28)
 hr_8= rep(8,28)


 hr = cbind(hr_0, hr_1, hr_2, hr_3, hr_4, hr_5, hr_6)

 hr=as.data.frame(hr)
 View(hr)


   Inc_Func = function (eir, start, Country, Region,
                               HR_sim, cov, min, max) {
   # First endectocide scenario
   ivm_parms = ivm_fun(IVM_start_times = start,
                        time_period = time_period,
                        hazard_profile = HR_sim,
                        ivm_coverage=cov,
                        ivm_min_age=min,
                        ivm_max_age = max)


   wh <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                     num_int = 2,
                                     het_brackets = 5,
                                     age = init_age,
                                     init_EIR = eir,
                                     country = Country,
                                     admin2 = Region,
                                     ttt = ivm_parms$ttt,
                                     eff_len = ivm_parms$eff_len,
                                     haz = ivm_parms$haz,
                                     ivm_cov_par = ivm_parms$ivm_cov_par,
                                     ivm_min_age = ivm_parms$ivm_min_age,
                                     ivm_max_age = ivm_parms$ivm_max_age,
                                     IVRM_start = ivm_parms$IVRM_start)


   res <- runfun(wh)

   res$time = res$t-3120

   return(res)
   }
View(hr)

##########
res0 = Inc_Func (eir=14, start=c(10000), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_0[1:28],cov=0.7, min=5, max=90)

res00 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                  HR_sim=hr$hr_0[1:28],cov=0.7, min=5, max=90)

res1 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_1[1:28],cov=0.7, min=5, max=90)

res2 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_2[1:28],cov=0.7, min=5, max=90)

res3 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_3[1:28],cov=0.7, min=5, max=90)

res4 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_4[1:28],cov=0.7, min=5, max=90)

res5 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_5[1:28],cov=0.7, min=5, max=90)

res6 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="=Senegal", Region="Fatick",
                 HR_sim=hr$hr_6[1:28],cov=0.7, min=5, max=90)

res7 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_7[1:28],cov=0.7, min=5, max=90)

res8 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                 HR_sim=hr$hr_8[1:28],cov=0.7, min=5, max=90)

########
   res0 = Inc_Func (eir=14, start=c(10000), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_0[1:28],cov=0.7, min=5, max=90)

   res00 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_0[1:28],cov=0.7, min=5, max=90)

   res1 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                              HR_sim=hr$hr_1[1:28],cov=0.7, min=5, max=90)

   res2 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                   HR_sim=hr$hr_2[1:28],cov=0.7, min=5, max=90)

   res3 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_3[1:28],cov=0.7, min=5, max=90)

   res4 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_4[1:28],cov=0.7, min=5, max=90)

   res5 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_5[1:28],cov=0.7, min=5, max=90)

   res6 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_6[1:28],cov=0.7, min=5, max=90)

   res7 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_7[1:28],cov=0.7, min=5, max=90)

   res8 = Inc_Func (eir=14,start=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim=hr$hr_8[1:28],cov=0.7, min=5, max=90)

   res00 = as.data.frame(res00)

   res0 = as.data.frame(res0)
   res1 = as.data.frame(res1)
   res2 = as.data.frame(res2)
   res3 = as.data.frame(res3)
   res4 = as.data.frame(res4)
   res5 = as.data.frame(res5)
   res6 = as.data.frame(res6)
   res7 = as.data.frame(res7)
   res8 = as.data.frame(res8)



   # Second endectocide scenario
   ivm_parms2 = ivm_fun(IVM_start_times = start2,
                        time_period = time_period,
                        hazard_profile = HR_sim2,
                        ivm_coverage=cov2,
                        ivm_min_age=min2,
                        ivm_max_age = max2)

   wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                     num_int = 2,
                                     het_brackets = 5,
                                     age = init_age,
                                     init_EIR = eir,
                                     country = Country,
                                     admin2 = Region,
                                     ttt = ivm_parms2$ttt,
                                     eff_len = ivm_parms2$eff_len,
                                     haz = ivm_parms2$haz,
                                     ivm_cov_par = ivm_parms2$ivm_cov_par,
                                     ivm_min_age = ivm_parms2$ivm_min_age,
                                     ivm_max_age = ivm_parms2$ivm_max_age,
                                     IVRM_start = ivm_parms2$IVRM_start)

   ggplot()+

     geom_line(data=res0, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="1"), col="grey", size=0.5) +
     geom_line(data=res00, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="2"), col='red',  size=0.5) +

     geom_line(data=res1, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="3"), col="blue", size=0.5) +
     geom_line(data=res2, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="4"), col="purple", size=0.5) +
     geom_line(data=res3, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="5"), col="orange", size=0.5) +
     geom_line(data=res4, aes(x=time/30, y=clin_inc0to80*1000*365, linetype="6"), col="brown", size=0.5) +

     xlab ("Time since start of intervention (Months)") +
     scale_x_continuous(expand = c(0,0), limits = c(-2,14), breaks=seq(-2,14,2)) +
     scale_y_continuous(expand = c(0,0), limits = c(0,1000)) +
     ylab("Incidence per 1,000 population")+
     scale_linetype_manual(values=c(3, 1, 1,1,1,1,1,1), labels = c("Basline malaria", "HR=0", "HR=1", "HR=2", "HR=3", "HR=4", "HR=5", "HR=6")) +

   guides(linetype = guide_legend(override.aes = list(color = c("grey", "red", "blue", "purple", "orange", "brown")))) +

     annotate("rect", xmin=0, xmax=1/28*3, ymin=0, ymax=900 ,  alpha=0.2, col="pink", fill="pink") +
     annotate("rect", xmin=1, xmax=1 + 1/28*3, ymin=0, ymax=900 , alpha=0.2, col="pink", fill="pink") +
     annotate("rect", xmin=2, xmax=2 + 1/28*3, ymin=0, ymax=900 , alpha=0.2, col="pink", fill="pink") +

     ggtitle ("Title") +
     theme(axis.text.x = element_text(size = 10),
           axis.text.y = element_text(size=10),
           axis.title.y = element_text(size =10),
           axis.title.x = element_text(size=10, vjust = -2),
           legend.title = element_blank(),legend.spacing.y = unit(1, 'mm'),
           axis.line = element_line(colour = "black", size = 0.6, linetype = "solid"),
           legend.position = c(0.7,0.9),
           #legend.box.background=element_rect(fill="white", size=0.5),
           legend.background= element_rect(fill="white"),
           legend.text = element_text(size=10),
           plot.title = element_text(size=10, hjust=0.5), legend.key.height =unit(0.1, 'cm'),
           #legend.box.margin = margin(0.1,0.1,0.1,0.1, unit="cm"), legend.margin = margin(-0.1,0,0,0, unit="cm"),
           panel.background = element_rect(fill="white"),
           panel.grid.major = element_line(colour = "white"),
           legend.key = element_rect(fill = alpha("white", 0.0)),
           legend.text.align = 0,
           plot.margin = margin(15,15,15,15),
           legend.key.width = unit(1, "cm"))


 ##########
   display.brewer.all()


 output_1 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=hr$hr_1[1:28],cov0=0.7, min0=5, max0=90,
                         HR_sim1=hr$hr_1[1:28], cov1=0.7, min1=5, max1=90,
                         HR_sim2=hr$hr_2[1:28], cov2=0.7, min2=5, max2=90)


 output_2 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=ivm_haz$IVM_300_3_AT_3[1:23],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$NTBC_1000_3_AT[1:23], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$NTBC_1000_3_AT[1:23], cov2=0.7, min2=0, max2=90)

 df = rbind(output_1, output_2)
 df = as.data.frame(df)
 colnames (df) = c("sum_inc0", "sum_inc1", "sum_inc2", "inc_red1", "inc_red2", "ave_prev0", "ave_prev1", "ave_prev2",  "prev_red1", "prev_red2")
 View(df)

  plot_1 = Inc_plot_Func_3 (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=hr$hr_1[1:28],cov0=0.7, min0=5, max0=90,
                           HR_sim1=hr$hr_3[1:28], cov1=0.7, min1=5, max1=90,
                           HR_sim2=hr$hr_8[1:28], cov2=0.7, min2=5, max2=90,
                           ylim1=c(0, 2200), ylim2=c(0,50), labels=c("IVM 3x300ug/kg (5-90yrs)", "IVM 1x400ug/kg(5-90yrs)", "Baseline malaria"),
                           Title= "Highly seasonal transmission")

 plot_2 = Inc_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT_3[1:23],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:23], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:23], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 2200), ylim2=c(0,50), labels=c("IVM 300ug/kg (5-90yrs)", "IVM 1x400ug/kg(5-90yrs)", "Baseline malaria"),
                           Title= "Perennial  transmission")


 plot_3 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                           HR_sim0=ivm_haz$IVM_300_3_AT_3[1:23],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:23], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:23], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("IVM 3x300ug/kg (5-90yrs)", "IVM 1x400ug/kg(5-90yrs)"),
                           Title= "Highly seasonal transmission")

 plot_4 = Inc_plot_Func_2 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                           HR_sim0=ivm_haz$IVM_300_3_AT_3[1:23],cov0=0.7, min0=5, max0=90,
                           HR_sim1=ivm_haz$NTBC_1000_3_AT[1:23], cov1=0.7, min1=5, max1=90,
                           HR_sim2=ivm_haz$NTBC_1000_3_AT[1:23], cov2=0.7, min2=0, max2=90,
                           ylim1=c(0, 900), ylim2=c(0,50), labels=c("IVM 3x300ug/kg (5-90yrs)", "IVM 1x400ug/kg(5-90yrs)"),
                           Title= "Perennial  transmission")


 plot_5 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                            HR_sim0=ivm_haz$IVM_300_3_AT_3[1:23],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_3_AT[1:23], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:23], cov2=0.7, min2=0, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg (5-90yrs)", "IVM 1x400ug/kg(5-90yrs)", "Baseline malaria"),
                            Title= "Highly seasonal transmission")


 plot_6 = Prev_plot_Func_3 (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                            HR_sim0=ivm_haz$IVM_300_3_AT_3[1:23],cov0=0.7, min0=5, max0=90,
                            HR_sim1=ivm_haz$NTBC_1000_3_AT[1:23], cov1=0.7, min1=5, max1=90,
                            HR_sim2=ivm_haz$NTBC_1000_3_AT[1:23], cov2=0.7, min2=0, max2=90,
                            ylim1=c(0, 2500), ylim2=c(0,50), labels=c("IVM 3x300ug/kg (5-90yrs)", "IVM 1x400ug/kg(5-90yrs)", "Baseline malaria"),
                            Title= "Perrenial transmission")





 ##############################


#Create dataframe of sinulation outputs
df = rbind(output_1, output_2)
df = as.data.frame(df)
colnames (df) = c("sum_inc0", "sum_inc1", "sum_inc2", "inc_red1", "inc_red2", "ave_prev0", "ave_prev1", "ave_prev2",  "prev_red1", "prev_red2")
View(df)


################
#Surplus code below


ivm_parms2 = ivm_fun(IVM_start_times = c(2250, 2280, 2310),
                     time_period = time_period,
                     hazard_profile = IVM_300[1:28],
                     ivm_coverage=cov2,
                     ivm_min_age=min2,
                     ivm_max_age = max2)

wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = eir,
                                  country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = Region,
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms2$ivm_cov_par,
                                  ivm_min_age = ivm_parms2$ivm_min_age,
                                  ivm_max_age = ivm_parms2$ivm_max_age,
                                  IVRM_start = ivm_parms2$IVRM_start)

# Running each of the models wh0 (no ivermectin), wh1 (ivermectin, HR = 2 for 10 days) and
# wh2 (ivermectin, HR = 2 for 28 days)
res0 <- runfun(wh0)
res1 <- runfun(wh1)
res2 <- runfun(wh2)

# Plotting the results
cols <- c("grey40", "deeppink2", "deepskyblue3")
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))

# Clinical Incidence
plot=plot(res0$t/365, res0$clin_inc0to80*1000, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 2, col = cols[1], xlim = c(0, 365), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$clin_inc0to80*1000, lwd = 2, col = cols[2])
lines(res2$t/365, res2$clin_inc0to80*1000, lwd = 2, col = cols[3])
#arrows(c(2250, 2280, 2310)/365, -50, c(2250, 2280, 2310)/365, 10, length = 0.15, lwd = 1, col = "goldenrod2")


#plot(res0$t/365, res0$slide_prev0to80*100, type="l", ylab="Slide prevalence (%)", xlab="Year", lwd=3, col=cols[1],
   #  xlim = c(5.7, 8), ylim = c(0, 70), las=1)
#lines(res1$t/365, res1$slide_prev0to80*100, lwd=3, col=cols[2])
#lines(res2$t/365, res2$slide_prev0to80*100, lwd=3, col=cols[3])
#arrows(c(2250, 2280, 2310)/365, -50, c(2250, 2280, 2310)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")

legend("topright", c("No endectocide", "IVM endectocide", "NTBC endectocide"),
       col = cols, lwd=3, bty="n", cex=0.8)
return (plot)

}

#########################################

ivm_parms0 = ivm_fun(IVM_start_times = 10000,
                     time_period = time_period,
                     hazard_profile = ivm_haz$IVM_300_3_AT[1:23],
                     ivm_coverage=0.8,
                     ivm_min_age=5,
                     ivm_max_age = 90)

wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = 14,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms0$ttt,
                                  eff_len = ivm_parms0$eff_len,
                                  haz = ivm_parms0$haz,
                                  ivm_cov_par = ivm_parms0$ivm_cov_par,
                                  ivm_min_age = ivm_parms0$ivm_min_age,
                                  ivm_max_age = ivm_parms0$ivm_max_age,
                                  IVRM_start = ivm_parms0$IVRM_start)



ivm_parms1 = ivm_fun(IVM_start_times = c(3120, 3150, 3180),
                     time_period = time_period,
                     hazard_profile = sim_ave_HR$IVM_300_3_AT[1:28],
                     ivm_coverage=0.8,
                     ivm_min_age=5,
                     ivm_max_age = 90)
wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = 14,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start)



ivm_parms2 = ivm_fun(IVM_start_times = c(3120, 3150, 3180),
                     time_period = time_period,
                     hazard_profile = sim_ave_HR$IVM_300_3_AT[1:28],
                     ivm_coverage=0.8,
                     ivm_min_age=0,
                     ivm_max_age = 90)

wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = 14,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms2$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms2$ivm_cov_par,
                                  ivm_min_age = ivm_parms2$ivm_min_age,
                                  ivm_max_age = ivm_parms2$ivm_max_age,
                                  IVRM_start = ivm_parms2$IVRM_start)

# Running each of the models wh0 (no ivermectin), wh1 (ivermectin, HR = 2 for 10 days) and
# wh2 (ivermectin, HR = 2 for 28 days)
res0 <- runfun(wh0)
res1 <- runfun(wh1)
res2 <- runfun(wh2)



sum_inc0=sum(res0$clin_inc0to80*1000*365)
sum_inc1=sum(res1$clin_inc0to80*1000*365)
sum_inc2=sum(res2$clin_inc0to80*1000*365)

inc_red1= sum_inc0 - sum_inc1
inc_red2= sum_inc0 - sum_inc2

sum_prev0=sum(res0$slide_prev0to80*100)
sum_prev1=sum(res1$slide_prev0to80*100)
sum_prev2=sum(res2$slide_prev0to80*100)

prev_red1= sum_prev0 - sum_prev1
prev_red2= sum_prev0 - sum_prev2

outcome_params <- c(sum_inc0, sum_inc1, sum_inc2, inc_red1, inc_red2, sum_prev0, sum_prev1,sum_prev2,  prev_red1, prev_red2)



# Plotting the results
cols <- c("grey40", "deeppink2", "deepskyblue3")
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))

# Clinical Incidence
plot(res0$t/365, res0$clin_inc0to80*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 2, col = cols[1], xlim = c(7.5, 10), ylim = c(0, 4000), las = 1)
lines(res1$t/365, res1$clin_inc0to80*1000*365, lwd = 2, col = cols[2])
lines(res2$t/365, res2$clin_inc0to80*1000*365, lwd = 2, col = cols[3])
arrows(c(3120, 3150, 3180)/365, -50, c(3120, 3150, 3180)/365, 10, length = 0.15, lwd = 1, col = "goldenrod2")


# Slide Prevalence
plot(res0$t/365, res0$slide_prev0to80*100, type="l", ylab="Slide prevalence (%)", xlab="Year", lwd=3, col=cols[1],
     xlim = c(7.5,10), ylim = c(0,70), las=1)
lines(res1$t/365, res1$slide_prev0to80*100, lwd=3, col=cols[2])
lines(res2$t/365, res2$slide_prev0to80*100, lwd=3, col=cols[3])
arrows(c(3120, 3150, 3180)/365, -50, c(3120, 3150, 3180)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")

legend("topright", c("No endectocide", "IVM endectocide", "NTBC endectocide"),
       col = cols, lwd=3, bty="n", cex=0.8)

}


# Testing Plotting All Outputs
par(mfrow = c(1, 3))

# Clinical Incidence 0-5 years, 0 to 59 months, 0-80 years
plot(res0$t/365, res0$clin_inc0to5*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$clin_inc0to5*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$clin_inc0to5*1000*365, lwd = 3, col = cols[3])

plot(res0$t/365, res0$clin_inc0to59*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$clin_inc0to59*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$clin_inc0to59*1000*365, lwd = 3, col = cols[3])


# Slide Prevalence 0-5 years, 0 to 59 months, 0-80 years
plot(res0$t/365, res0$slide_prev0to5, type = "l", ylab = "Slide prevalence (%)", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1), las = 1)
lines(res1$t/365, res1$slide_prev0to5, lwd = 3, col = cols[2])
lines(res2$t/365, res2$slide_prev0to5, lwd = 3, col = cols[3])

plot(res0$t/365, res0$slide_prev0to59, type = "l", ylab = "Slide prevalence (%)", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1), las = 1)
lines(res1$t/365, res1$slide_prev0to59, lwd = 3, col = cols[2])
lines(res2$t/365, res2$slide_prev0to59, lwd = 3, col = cols[3])

plot(res0$t/365, res0$slide_prev0to80, type = "l", ylab = "Slide prevalence (%)", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1), las = 1)
lines(res1$t/365, res1$slide_prev0to80, lwd = 3, col = cols[2])
lines(res2$t/365, res2$slide_prev0to80, lwd = 3, col = cols[3])

# Deaths 0-5 years, 0 to 59 months, 0-80 years
plot(res0$t/365, res0$deaths_inc0to5*1000*365, type = "l", ylab = "Deaths", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$deaths_inc0to5*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$deaths_inc0to5*1000*365, lwd = 3, col = cols[3])

plot(res0$t/365, res0$deaths_inc0to59*1000*365, type = "l", ylab = "Deaths", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$deaths_inc0to59*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$deaths_inc0to59*1000*365, lwd = 3, col = cols[3])

plot(res0$t/365, res0$deaths_inc0to80*1000*365, type = "l", ylab = "Deaths", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$deaths_inc0to80*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$deaths_inc0to80*1000*365, lwd = 3, col = cols[3])





