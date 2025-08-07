require(tidyverse)
# i can't easily find the 1x400 HR numbers so i've made up something close-ish
t <- 1:30
hr <- c(1, 9, 1+8*exp(-0.2*(t-2)))
plot(t, hr[1:30], type = "o", pch=16)
abline(h=1, lty=2)


# function to generate hazard curves based on when each individual recieved IVM in the MDA
#NC: I think npop is equivalent to number of distributions.
make_haz_curves <- function(ivm_mda_duration = ivm_mda_duration,
                            ivm_effect_duration = 23, #for IVM HRs
                            npop = npop){


  all_haz <- NULL
  tseq <- 1:(ivm_mda_duration + ivm_effect_duration) #time goes from duration of MDA + max time of effect --> have enough time steps

  for(i in 1:npop){

    # randomly select the day on which individual i took ivermectin
    ivm_start_time <- round(runif(1, 1, ivm_mda_duration),0) #randomly draw a start time within the ivm_mda_duration

    # make an empty dataframe for that individuals hr data
    haz_curve <- data.frame(ind = i,
                            tt = tseq,
                            haz = 1) #at the beginning, no killing effect.
    # fill in the haz curve at the timepoint they took ivm
    haz_curve$haz[ivm_start_time:(ivm_start_time + ivm_effect_duration-1)] = hr[1:ivm_effect_duration]

    all_haz <- bind_rows(all_haz, haz_curve)
  }

  return(all_haz)
}

#### test the function works ####

pop_haz = make_haz_curves(ivm_mda_duration = 10,
                npop = 100)

ggplot(pop_haz, aes(x = tt, y = haz, color = as.factor(ind))) +
  geom_line() +
  theme(legend.position = "none")

# then take the average hazard at each day
mean_daily_hazard <- pop_haz %>%
  group_by(time_since_mda_started = as.factor(tt)) %>%
  summarise(mean_daily_hazard = mean(haz))


ggplot(mean_daily_hazard, aes(x = time_since_mda_started, y = mean_daily_hazard)) +
  geom_col()

#### comparison for 11 days vs 31 days ####


pop_haz_11 = make_haz_curves(ivm_mda_duration = 11,
                          npop = 10000)

pop_haz_31 = make_haz_curves(ivm_mda_duration = 31,
                             npop = 10000)


mean_daily_hazard_11 <- pop_haz_11 %>%
  group_by(time_since_mda_started = as.factor(tt)) %>%
  summarise(mean_daily_hazard = mean(haz)) %>%
  mutate(mda_duration = 11)

mean_daily_hazard_31 <- pop_haz_31 %>%
  group_by(time_since_mda_started = as.factor(tt)) %>%
  summarise(mean_daily_hazard = mean(haz)) %>%
  mutate(mda_duration = 31)


combined_data <- bind_rows(mean_daily_hazard_11, mean_daily_hazard_31)

ggplot(combined_data, aes(x = as.numeric(time_since_mda_started),
                          y = mean_daily_hazard,
                          color = as.factor(mda_duration))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1) +
  ylim(c(0, 7)) +
  theme_bw()

