

setwd("~/OneDrive - LSTM/ivRmectin")
# Loading the ivRmectin package
devtools::load_all()


# Create a vector of age categories for the model
init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 3

#low - 2, moderate - 15, high - 120


# Provide the length of time (in days) that you want to run the model for
time_period <- 365 *10 # run model for 10 years

# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")
library(RColorBrewer)

# Load ivermectin hazardd
ivm_haz_1 <- read.table("data/ivermectin_hazards.txt", header=TRUE)
ivm_haz_2 <- read.csv("data/sim_ave_HR.csv", header=TRUE)

ivm_haz_2 <- read.csv("/Users/anna.trett/OneDrive - LSTM/CHIC 599/NTBC_Vs_IVM_Summer2020/PKProfiles/01-Data/sim_ave_HR.csv", header=TRUE)

ivm_haz_1 = ivm_haz_1[1:23,]
ivm_haz_2 = ivm_haz_2[1:23, ]
ivm_haz=cbind(ivm_haz_1, ivm_haz_2)
View(ivm_haz)

library("dplyr")
ivm_haz = ivm_haz %>%
  select (-"Day")

colnames(ivm_haz) = c("Day", "IVM_400_1", "IVM_300_3", "NTBC_1000_3", "NTBC_500_3", "NTBC_1000_1")



# Running the ivm_fun function to generate the extra endectocide specific parameters that
# you have to pass to the model
ivm_parms = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                     time_period = time_period, # time period for the model to run over
                     hazard_profile = ivm_haz$IVM_300_3[1:23], # dummy hazard profile - must be vector
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



#loop to obtain optimal timing of seasonal MDA

##Function - Optimisation

vector=vector(mode="numeric",length=3)
list=vector(mode="list",length=10)


for (i in 1:1215) {
  vector[1] <- 30*i
  vector[2] <- vector[1] + 30
  vector[3] <- vector[2] + 30
  list[[i]] <- vector
}

list[104]
sum_inc_list=vector(mode="list",length=1215)

for (i in 1:1215) {
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
                       hazard_profile = ivm_haz$IVM_400_1, # dummy hazard profile - must be vector
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

  res <- runfun(wh0)

  sum_inc=sum(res$clin_inc0to80*1000*365)

  sum_inc_list[i] <- sum_inc

}

sum_inc_list = exciting

administration = unlist(sum_inc_list)

View(administration)
min(administration)

optim(list), Sum_incidence_Func


#################


###############################################
#Function for generating clinical incidence and slide prevalence reductions for three scenarios
Params_Func = function (eir, start1, start2, Country, Region,
                        HR_sim0, cov0, min0, max0,
                        HR_sim1, cov1, min1, max1,
                        HR_sim2, cov2, min2, max2,
                        ylim1, ylim2, labels) {

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

  # Plotting the results


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

  return (outcome_params)

}

#Function which returns total incidence and prevalence statistics for comparison
Plot_Func = function (eir, start1, start2, Country, Region,
                        HR_sim0, cov0, min0, max0,
                        HR_sim1, cov1, min1, max1,
                        HR_sim2, cov2, min2, max2,
                        ylim1, ylim2, labels) {

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

  # Plotting the results
  cols <- c("grey40", "deeppink2", "deepskyblue3")
  par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))

  res0$time = res0$t-3120
  res1$time = res1$t-3120
  res2$time = res2$t-3120


  ggplot()+

    # geom_line(data=res0, aes(x=time/120, y=clin_inc0to80*1000*365, linetype="3"), color="darkgrey", size=0.5) +

    geom_line(data=res1, aes(x=time/28, y=clin_inc0to80*1000*365, linetype="1"), color="#689ed4", size=0.5) +
    geom_line(data=res2, aes(x=time/28, y=clin_inc0to80*1000*365, linetype="2"), color="#c07142", size=0.5) +

    xlab ("Time since start of intervention (Months)") +
    scale_x_continuous(expand = c(0,0), limits = c(-2,6), breaks=seq(-2,6,2)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,700)) +
    ylab("Incidence per 1,000 population")+
    scale_linetype_manual(values=c(1, 1), labels = c("IVM 300ug/kg", "NTBC 1mg/kg")) +
    guides(linetype = guide_legend(override.aes = list(color = c("#689ed4", "#c07142")))) +

    annotate("rect", xmin=0, xmax=1/28*3, ymin=0, ymax=700 ,  col="pink") +
    annotate("rect", xmin=1, xmax=1 + 1/28*3, ymin=0, ymax=700 , col="pink") +
    annotate("rect", xmin=2, xmax=2 + 1/28*3, ymin=0, ymax=700 , col="pink") +

    ggtitle ("Highly seasonal transmission") +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10),
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



    ggplot()+

    # geom_line(data=res0, aes(x=time/120, y=clin_inc0to80*1000*365, linetype="3"), color="darkgrey", size=0.5) +

    geom_line(data=res1, aes(x=time/28, y=slide_prev0to80*100, linetype="1"), color="#689ed4", size=0.5) +
    geom_line(data=res2, aes(x=time/28, y=slide_prev0to80*100, linetype="2"), color="#c07142", size=0.5) +

    xlab ("Time since start of intervention (Months)") +
    scale_x_continuous(expand = c(0,0), limits = c(-2,6), breaks=seq(-2,6,2)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,40)) +
    ylab("Incidence per 1,000 population")+
    scale_linetype_manual(values=c(1, 1), labels = c("IVM 300ug/kg", "NTBC 1mg/kg")) +
    guides(linetype = guide_legend(override.aes = list(color = c("#689ed4", "#c07142")))) +

    annotate("rect", xmin=0, xmax=1/28*3, ymin=0, ymax=40 ,  col="pink") +
    annotate("rect", xmin=1, xmax=1 + 1/28*3, ymin=0, ymax=40 , col="pink") +
    annotate("rect", xmin=2, xmax=2 + 1/28*3, ymin=0, ymax=40 , col="pink") +

    ggtitle ("Highly seasonal transmission") +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10),
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
  ########################################


  ggplot()+

   # geom_line(data=res0, aes(x=time/120, y=clin_inc0to80*1000*365, linetype="3"), color="darkgrey", size=0.5) +

    geom_line(data=res1, aes(x=time/28, y=clin_inc0to80*1000*365, linetype="1"), color="#689ed4", size=0.5) +
    geom_line(data=res2, aes(x=time/28, y=clin_inc0to80*1000*365, linetype="2"), color="#c07142", size=0.5) +

    xlab ("Time since start of intervention (Months)") +
    scale_x_continuous(expand = c(0,0), limits = c(-2,6), breaks=seq(-2,6,2)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,700)) +
    ylab("Incidence per 1,000 population")+
    scale_linetype_manual(values=c(1, 1), labels = c("IVM 300ug/kg", "NTBC 1mg/kg")) +
    guides(linetype = guide_legend(override.aes = list(color = c("#689ed4", "#c07142")))) +

    annotate("rect", xmin=0, xmax=1/28*3, ymin=0, ymax=700 ,  col="pink") +
    annotate("rect", xmin=1, xmax=1 + 1/28*3, ymin=0, ymax=700 , col="pink") +
    annotate("rect", xmin=2, xmax=2 + 1/28*3, ymin=0, ymax=700 , col="pink") +

    ggtitle ("Highly seasonal transmission") +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size=10),
          axis.title.y = element_text(size=10),
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

  scale_color_brewer()

library(RColorBrewer)

  # Clinical Incidence
  plot(res0$t/365, res0$clin_inc0to80*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
       lwd = 1, col = cols[1], xlim = c(7.5, 10), ylim = ylim1, las = 1)
  lines(res1$t/365, res1$clin_inc0to80*1000*365, lwd = 1, col = cols[2])
  lines(res2$t/365, res2$clin_inc0to80*1000*365, lwd = 1, col = cols[3])
  arrows(c(3120, 3150, 3180)/365, -200, c(3120, 3150, 3180)/365, 1, length = 0.1, lwd = 1, col = "goldenrod2")

  legend("topright", labels,
         col = cols, lwd=2, bty="n", cex=0.8)

  plot(res0$t/365, res0$slide_prev0to80*100, type = "l", ylab ="Slide prevalence (%)", xlab = "Year",
       lwd = 1, col = cols[1], xlim = c(7.5, 10), ylim = ylim2, las = 1)
  lines(res1$t/365, res1$slide_prev0to80*100, lwd = 1, col = cols[2])
  lines(res2$t/365, res2$slide_prev0to80*100, lwd = 1, col = cols[3])
  arrows(c(3120, 3150, 3180)/365, -50, c(3120, 3150, 3180)/365, 1, length = 0.1, lwd = 1, col = "goldenrod2")

  legend("topright", labels,
         col = cols, lwd=2, bty="n", cex=0.8)

  return (plot)

}


###############
#graphical output

#Comparison of age coverage between IVM maximum dose
#Adult incidence and prevalence outcomes
#highly seasonal
output_1 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_300_3[1:23], cov2=0.7, min2=0, max2=90,
                         ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))



plot_1 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                         HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                         HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                         HR_sim2=ivm_haz$IVM_300_3[1:23], cov2=0.7, min2=0, max2=90,
                         ylim1=c(0, 2500), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))


output_2= Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                        HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                        HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                        HR_sim2=ivm_haz$IVM_300_3[1:23], cov2=0.7, min2=0, max2=90,
                        ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))



plot_2 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                    HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                    HR_sim2=ivm_haz$IVM_300_3[1:23], cov2=0.7, min2=0, max2=90,
                    ylim1=c(0, 1000), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))


###################################3
output_1 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                        HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                        HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                        HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=5, max2=90,
                        ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs"))



plot_1 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                    HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                    HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                    HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=5, max2=90,
                    ylim1=c(0, 2500), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))


output_2= Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                       HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                       HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                       HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=5, max2=90,
                       ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))



plot_2 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                    HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                    HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=5, max2=90,
                    ylim1=c(0, 1000), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))


####################################


output_1 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                        HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                        HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                        HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                        ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-5yrs)"))



plot_1 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                    HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                    HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                    HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                    ylim1=c(0, 2500), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-5yrs)"))


output_2= Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                       HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                       HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                       HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                       ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-5yrs)"))



plot_2 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                    HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                    HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                    ylim1=c(0, 1000), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-5yrs)"))




####################################


output_1 = Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                        HR_sim0=ivm_haz$NTBC_1000_1[1:23],cov0=0.7, min0=5, max0=90,
                        HR_sim1=ivm_haz$NTBC_1000_1[1:23], cov1=0.7, min1=5, max1=90,
                        HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                        ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "NTBC (3 doses)", "NTBC (1 dose)"))



plot_1 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Senegal", Region="Fatick",
                    HR_sim0=ivm_haz$NTBC_1000_1[1:23],cov0=0.7, min0=5, max0=90,
                    HR_sim1=ivm_haz$NTBC_1000_1[1:23], cov1=0.7, min1=5, max1=90,
                    HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                    ylim1=c(0, 2500), ylim2=c(0,70), labels=c("No endectocide", "NTBC (3 doses)", "NTBC (1 dose)"))


output_2= Params_Func (eir=14, start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                       HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                       HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                       HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                       ylim1=c(0, 2500), ylim2= c(0,60), labels=c("No endectocide", "NTBC (3 doses)", "NTBC (1 dose)"))



plot_2 = Plot_Func (eir=14,   start1=c(3120, 3150, 3180), start2=c(3120, 3150, 3180), Country="Democratic Republic of the Congo", Region="Equateur",
                    HR_sim0=ivm_haz$IVM_300_3[1:23],cov0=0.7, min0=5, max0=90,
                    HR_sim1=ivm_haz$IVM_300_3[1:23], cov1=0.7, min1=5, max1=90,
                    HR_sim2=ivm_haz$NTBC_1000_3[1:23], cov2=0.7, min2=0, max2=90,
                    ylim1=c(0, 1000), ylim2=c(0,70), labels=c("No endectocide", "NTBC (3 doses)", "NTBC (1 dose)"))










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
                     hazard_profile = ivm_haz$IVM_300[1:23],
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
                     hazard_profile = sim_ave_HR$IVM_300[1:28],
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
                     hazard_profile = sim_ave_HR$IVM_300[1:28],
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





