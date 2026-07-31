#find best-fitting values for endec_mu and wane for each dosage and setting

#kwale 1x400

mda_start <- (30*3)
IVM_begin1 <- (365*12)+mda_start

mda_int <- 30 #every 30 days

IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

my_sim_kwale_out_1_400 <- readRDS("analysis/IVER_ONE/output/out_kwale_1_400.rds")
my_sim_kwale_out_1_400_exp_decay <- readRDS("analysis/IVER_ONE/output/out_kwale_1_400_exp_decay.rds")
my_sim_kwale_out_1_800 <- readRDS("analysis/IVER_ONE/output/out_kwale_1_800.rds")
my_sim_kwale_out_1_800_exp_decay <- readRDS("analysis/IVER_ONE/output/out_kwale_1_800_exp_decay.rds")
my_sim_kwale_out_3_300 <- readRDS("analysis/IVER_ONE/output/out_kwale_3_300.rds")
my_sim_kwale_out_3_300_exp_decay <- readRDS("analysis/IVER_ONE/output/out_kwale_3_300_exp_decay.rds")
my_sim_kwale_out_3_600 <- readRDS("analysis/IVER_ONE/output/out_kwale_3_600.rds")
my_sim_kwale_out_3_600_exp_decay <- readRDS("analysis/IVER_ONE/output/out_kwale_3_600_exp_decay.rds")

#busia
my_sim_busia_out_1_400 <- readRDS("analysis/IVER_ONE/output/out_busia_1_400.rds")
my_sim_busia_out_1_400_exp_decay <- readRDS("analysis/IVER_ONE/output/out_busia_1_400_exp_decay.rds")
my_sim_busia_out_1_800 <- readRDS("analysis/IVER_ONE/output/out_busia_1_800.rds")
my_sim_busia_out_1_800_exp_decay <- readRDS("analysis/IVER_ONE/output/out_busia_1_800_exp_decay.rds")
my_sim_busia_out_3_300 <- readRDS("analysis/IVER_ONE/output/out_busia_3_300.rds")
my_sim_busia_out_3_300_exp_decay <- readRDS("analysis/IVER_ONE/output/out_busia_3_300_exp_decay.rds")
my_sim_busia_out_3_600 <- readRDS("analysis/IVER_ONE/output/out_busia_3_600.rds")
my_sim_busia_out_3_600_exp_decay <- readRDS("analysis/IVER_ONE/output/out_busia_3_600_exp_decay.rds")

#homa bay
my_sim_homa_bay_out_1_400 <- readRDS("analysis/IVER_ONE/output/out_hb_1_400.rds")
my_sim_homa_bay_out_1_400_exp_decay <- readRDS("analysis/IVER_ONE/output/out_hb_1_400_exp_decay.rds")
my_sim_homa_bay_out_1_800 <- readRDS("analysis/IVER_ONE/output/out_hb_1_800.rds")
my_sim_homa_bay_out_1_800_exp_decay <- readRDS("analysis/IVER_ONE/output/out_hb_1_800_exp_decay.rds")
my_sim_homa_bay_out_3_300 <- readRDS("analysis/IVER_ONE/output/out_hb_3_300.rds")
my_sim_homa_bay_out_3_300_exp_decay <- readRDS("analysis/IVER_ONE/output/out_hb_3_300_exp_decay.rds")
my_sim_homa_bay_out_3_600 <- readRDS("analysis/IVER_ONE/output/out_hb_3_600.rds")
my_sim_homa_bay_out_3_600_exp_decay <- readRDS("analysis/IVER_ONE/output/out_hb_3_600_exp_decay.rds")

#migori

my_sim_migori_out_1_400 <- readRDS("analysis/IVER_ONE/output/out_migori_1_400.rds")
my_sim_migori_out_1_400_exp_decay <- readRDS("analysis/IVER_ONE/output/out_migori_1_400_exp_decay.rds")
my_sim_migori_out_1_800 <- readRDS("analysis/IVER_ONE/output/out_migori_1_800.rds")
my_sim_migori_out_1_800_exp_decay <- readRDS("analysis/IVER_ONE/output/out_migori_1_800_exp_decay.rds")
my_sim_migori_out_3_300 <- readRDS("analysis/IVER_ONE/output/out_migori_3_300.rds")
my_sim_migori_out_3_300_exp_decay <- readRDS("analysis/IVER_ONE/output/out_migori_3_300_exp_decay.rds")
my_sim_migori_out_3_600 <- readRDS("analysis/IVER_ONE/output/out_migori_3_600.rds")
my_sim_migori_out_3_600_exp_decay <- readRDS("analysis/IVER_ONE/output/out_migori_3_600_exp_decay.rds")

#siaya
my_sim_siaya_out_1_400 <- readRDS("analysis/IVER_ONE/output/out_siaya_1_400.rds")
my_sim_siaya_out_1_400_exp_decay <- readRDS("analysis/IVER_ONE/output/out_siaya_1_400_exp_decay.rds")
my_sim_siaya_out_1_800 <- readRDS("analysis/IVER_ONE/output/out_siaya_1_800.rds")
my_sim_siaya_out_1_800_exp_decay <- readRDS("analysis/IVER_ONE/output/out_siaya_1_800_exp_decay.rds")
my_sim_siaya_out_3_300 <- readRDS("analysis/IVER_ONE/output/out_siaya_3_300.rds")
my_sim_siaya_out_3_300_exp_decay <- readRDS("analysis/IVER_ONE/output/out_siaya_3_300_exp_decay.rds")
my_sim_siaya_out_3_600 <- readRDS("analysis/IVER_ONE/output/out_siaya_3_600.rds")
my_sim_siaya_out_3_600_exp_decay <- readRDS("analysis/IVER_ONE/output/out_siaya_3_600_exp_decay.rds")
####

endec_t_start_1_400 <- IVM_start[1]#
endec_t_end_1_400 <- IVM_start[3]+58 ##update

endec_t_start <- IVM_start[1]
endec_t_end <- IVM_start[3]+60 ##update

#filter both models to time when endec is active
times_kwale_out_1_400_exp_decay <- my_sim_kwale_out_1_400_exp_decay %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

times_kwale_odin_1_400 <- my_sim_kwale_out_1_400 %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_kwale_1_400 <- split(times_kwale_out_1_400_exp_decay, f = times_kwale_out_1_400_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_kwale_odin_1_400$Ivtot - list_kwale_1_400[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_kwale_1_400 <- malsim_odin[index,]

fit_kwale_1_400 <- my_sim_kwale_out_1_400_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_kwale_out_1_400, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_kwale_1_400, aes(x = t, y = Ivtot), col = "red")

#1_800
times_kwale_out_1_800_exp_decay <- my_sim_kwale_out_1_800_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_kwale_odin_1_800 <- my_sim_kwale_out_1_800 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_kwale_1_800 <- split(times_kwale_out_1_800_exp_decay, f = times_kwale_out_1_800_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_kwale_odin_1_800$Ivtot - list_kwale_1_800[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_kwale_1_800 <- malsim_odin[index,]

fit_kwale_1_800 <- my_sim_kwale_out_1_800_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_kwale_out_1_800, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_kwale_1_800, aes(x = t, y = Ivtot), col = "red")

#3-300
times_kwale_out_3_300_exp_decay <- my_sim_kwale_out_3_300_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_kwale_odin_3_300 <- my_sim_kwale_out_3_300 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_kwale_3_300 <- split(times_kwale_out_3_300_exp_decay, f = times_kwale_out_3_300_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_kwale_odin_3_300$Ivtot - list_kwale_3_300[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_kwale_3_300 <- malsim_odin[index,]

fit_kwale_3_300 <- my_sim_kwale_out_3_300_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_kwale_out_3_300, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_kwale_3_300, aes(x = t, y = Ivtot), col = "red")

#3-600
times_kwale_out_3_600_exp_decay <- my_sim_kwale_out_3_600_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_kwale_odin_3_600 <- my_sim_kwale_out_3_600 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_kwale_3_600 <- split(times_kwale_out_3_600_exp_decay, f = times_kwale_out_3_600_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_kwale_odin_3_600$Ivtot - list_kwale_3_600[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_kwale_3_600 <- malsim_odin[index,]

fit_kwale_3_600 <- my_sim_kwale_out_3_600_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_kwale_out_3_600, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_kwale_3_600, aes(x = t, y = Ivtot), col = "red")

#sort output
site <- "kwale"
dose <- c("1x400", "1x800", "3x300", "3x600")
endec_mu_vec <- c(best_fit_kwale_1_400$endec_mu, best_fit_kwale_1_800$endec_mu,
                  best_fit_kwale_3_300$endec_mu,
                  best_fit_kwale_3_600$endec_mu)

wane_vec <- c(best_fit_kwale_1_400$wane, best_fit_kwale_1_800$wane,
              best_fit_kwale_3_300$wane,
              best_fit_kwale_3_600$wane)

df_kwale <- data.frame(site, dose, endec_mu_vec, wane_vec)

#busia

#1x400
times_busia_out_1_400_exp_decay <- my_sim_busia_out_1_400_exp_decay %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

times_busia_odin_1_400 <- my_sim_busia_out_1_400 %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_busia_1_400 <- split(times_busia_out_1_400_exp_decay, f = times_busia_out_1_400_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_busia_odin_1_400$Ivtot - list_busia_1_400[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_busia_1_400 <- malsim_odin[index,]

fit_busia_1_400 <- my_sim_busia_out_1_400_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_busia_out_1_400, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_busia_1_400, aes(x = t, y = Ivtot), col = "red")

#1x800
times_busia_out_1_800_exp_decay <- my_sim_busia_out_1_800_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_busia_odin_1_800 <- my_sim_busia_out_1_800 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_busia_1_800 <- split(times_busia_out_1_800_exp_decay, f = times_busia_out_1_800_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_busia_odin_1_800$Ivtot - list_busia_1_800[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_busia_1_800 <- malsim_odin[index,]

fit_busia_1_800 <- my_sim_busia_out_1_800_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_busia_out_1_800, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_busia_1_800, aes(x = t, y = Ivtot), col = "red")

#3x300
times_busia_out_3_300_exp_decay <- my_sim_busia_out_3_300_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_busia_odin_3_300 <- my_sim_busia_out_3_300 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_busia_3_300 <- split(times_busia_out_3_300_exp_decay, f = times_busia_out_3_300_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_busia_odin_3_300$Ivtot - list_busia_3_300[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_busia_3_300 <- malsim_odin[index,]

fit_busia_3_300 <- my_sim_busia_out_3_300_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_busia_out_3_300, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_busia_3_300, aes(x = t, y = Ivtot), col = "red")

#3x600
times_busia_out_3_600_exp_decay <- my_sim_busia_out_3_600_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_busia_odin_3_600 <- my_sim_busia_out_3_600 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_busia_3_600 <- split(times_busia_out_3_600_exp_decay, f = times_busia_out_3_600_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_busia_odin_3_600$Ivtot - list_busia_3_600[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_busia_3_600 <- malsim_odin[index,]

fit_busia_3_600 <- my_sim_busia_out_3_600_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_busia_out_3_600, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_busia_3_600, aes(x = t, y = Ivtot), col = "red")

site <- "busia"
dose <- c("1x400", "1x800", "3x300", "3x600")
endec_mu_vec <- c(best_fit_busia_1_400$endec_mu, best_fit_busia_1_800$endec_mu,
                  best_fit_busia_3_300$endec_mu,
                  best_fit_busia_3_600$endec_mu)

wane_vec <- c(best_fit_busia_1_400$wane, best_fit_busia_1_800$wane,
              best_fit_busia_3_300$wane,
              best_fit_busia_3_600$wane)

df_busia <- data.frame(site, dose, endec_mu_vec, wane_vec)

#homa bay
#1x400
times_hb_out_1_400_exp_decay <- my_sim_homa_bay_out_1_400_exp_decay %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

times_hb_odin_1_400 <- my_sim_homa_bay_out_1_400 %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_hb_1_400 <- split(times_hb_out_1_400_exp_decay, f = times_hb_out_1_400_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_hb_odin_1_400$Ivtot - list_hb_1_400[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_hb_1_400 <- malsim_odin[index,]

fit_hb_1_400 <- my_sim_homa_bay_out_1_400_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_homa_bay_out_1_400, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_hb_1_400, aes(x = t, y = Ivtot), col = "red")

#1x800
times_hb_out_1_800_exp_decay <- my_sim_homa_bay_out_1_800_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_hb_odin_1_800 <- my_sim_homa_bay_out_1_800 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_hb_1_800 <- split(times_hb_out_1_800_exp_decay, f = times_hb_out_1_800_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_hb_odin_1_800$Ivtot - list_hb_1_800[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_hb_1_800 <- malsim_odin[index,]

fit_hb_1_800 <- my_sim_homa_bay_out_1_800_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_homa_bay_out_1_800, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_hb_1_800, aes(x = t, y = Ivtot), col = "red")

#3x300
times_hb_out_3_300_exp_decay <- my_sim_homa_bay_out_3_300_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_hb_odin_3_300 <- my_sim_homa_bay_out_3_300 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_hb_3_300 <- split(times_hb_out_3_300_exp_decay, f = times_hb_out_3_300_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_hb_odin_3_300$Ivtot - list_hb_3_300[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_hb_3_300 <- malsim_odin[index,]

fit_hb_3_300 <- my_sim_homa_bay_out_3_300_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_homa_bay_out_3_300, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_hb_3_300, aes(x = t, y = Ivtot), col = "red")

#3x600
times_hb_out_3_600_exp_decay <- my_sim_homa_bay_out_3_600_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_hb_odin_3_600 <- my_sim_homa_bay_out_3_600 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_hb_3_600 <- split(times_hb_out_3_600_exp_decay, f = times_hb_out_3_600_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_hb_odin_3_600$Ivtot - list_hb_3_600[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_hb_3_600 <- malsim_odin[index,]

fit_hb_3_600 <- my_sim_homa_bay_out_3_600_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_homa_bay_out_3_600, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_hb_3_600, aes(x = t, y = Ivtot), col = "red")

site <- "homa bay"
dose <- c("1x400", "1x800", "3x300", "3x600")
endec_mu_vec <- c(best_fit_hb_1_400$endec_mu, best_fit_hb_1_800$endec_mu,
                  best_fit_hb_3_300$endec_mu,
                  best_fit_hb_3_600$endec_mu)

wane_vec <- c(best_fit_hb_1_400$wane, best_fit_hb_1_800$wane,
              best_fit_hb_3_300$wane,
              best_fit_hb_3_600$wane)

df_hb <- data.frame(site, dose, endec_mu_vec, wane_vec)

#MIGORI
#1x400
times_migori_out_1_400_exp_decay <- my_sim_migori_out_1_400_exp_decay %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

times_migori_odin_1_400 <- my_sim_migori_out_1_400 %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_migori_1_400 <- split(times_migori_out_1_400_exp_decay, f = times_migori_out_1_400_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_migori_odin_1_400$Ivtot - list_migori_1_400[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_migori_1_400 <- malsim_odin[index,]

fit_migori_1_400 <- my_sim_migori_out_1_400_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_migori_out_1_400, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_migori_1_400, aes(x = t, y = Ivtot), col = "red")

#1x800
times_migori_out_1_800_exp_decay <- my_sim_migori_out_1_800_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_migori_odin_1_800 <- my_sim_migori_out_1_800 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_migori_1_800 <- split(times_migori_out_1_800_exp_decay, f = times_migori_out_1_800_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_migori_odin_1_800$Ivtot - list_migori_1_800[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_migori_1_800 <- malsim_odin[index,]

fit_migori_1_800 <- my_sim_migori_out_1_800_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_migori_out_1_800, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_migori_1_800, aes(x = t, y = Ivtot), col = "red")

#3x300
times_migori_out_3_300_exp_decay <- my_sim_migori_out_3_300_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_migori_odin_3_300 <- my_sim_migori_out_3_300 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_migori_3_300 <- split(times_migori_out_3_300_exp_decay, f = times_migori_out_3_300_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_migori_odin_3_300$Ivtot - list_migori_3_300[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_migori_3_300 <- malsim_odin[index,]

fit_migori_3_300 <- my_sim_migori_out_3_300_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_migori_out_3_300, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_migori_3_300, aes(x = t, y = Ivtot), col = "red")

#3x600
times_migori_out_3_600_exp_decay <- my_sim_migori_out_3_600_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_migori_odin_3_600 <- my_sim_migori_out_3_600 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_migori_3_600 <- split(times_migori_out_3_600_exp_decay, f = times_migori_out_3_600_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_migori_odin_3_600$Ivtot - list_migori_3_600[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_migori_3_600 <- malsim_odin[index,]

fit_migori_3_600 <- my_sim_migori_out_3_600_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_migori_out_3_600, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_migori_3_600, aes(x = t, y = Ivtot), col = "red")

site <- "migori"
dose <- c("1x400", "1x800", "3x300", "3x600")
endec_mu_vec <- c(best_fit_migori_1_400$endec_mu, best_fit_migori_1_800$endec_mu,
                  best_fit_migori_3_300$endec_mu,
                  best_fit_migori_3_600$endec_mu)

wane_vec <- c(best_fit_migori_1_400$wane, best_fit_migori_1_800$wane,
              best_fit_migori_3_300$wane,
              best_fit_migori_3_600$wane)

df_migori <- data.frame(site, dose, endec_mu_vec, wane_vec)

#SIAYA
#1x400

times_siaya_out_1_400_exp_decay <- my_sim_siaya_out_1_400_exp_decay %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

times_siaya_odin_1_400 <- my_sim_siaya_out_1_400 %>%
  filter(between(t, endec_t_start_1_400, endec_t_end_1_400))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_siaya_1_400 <- split(times_siaya_out_1_400_exp_decay, f = times_siaya_out_1_400_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_siaya_odin_1_400$Ivtot - list_siaya_1_400[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_siaya_1_400 <- malsim_odin[index,]

fit_siaya_1_400 <- my_sim_siaya_out_1_400_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_siaya_out_1_400, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_siaya_1_400, aes(x = t, y = Ivtot), col = "red")

#1x800
times_siaya_out_1_800_exp_decay <- my_sim_siaya_out_1_800_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_siaya_odin_1_800 <- my_sim_siaya_out_1_800 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_siaya_1_800 <- split(times_siaya_out_1_800_exp_decay, f = times_siaya_out_1_800_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_siaya_odin_1_800$Ivtot - list_siaya_1_800[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_siaya_1_800 <- malsim_odin[index,]

fit_siaya_1_800 <- my_sim_siaya_out_1_800_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_siaya_out_1_800, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_siaya_1_800, aes(x = t, y = Ivtot), col = "red")

#3x300
times_siaya_out_3_300_exp_decay <- my_sim_siaya_out_3_300_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_siaya_odin_3_300 <- my_sim_siaya_out_3_300 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_siaya_3_300 <- split(times_siaya_out_3_300_exp_decay, f = times_siaya_out_3_300_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_siaya_odin_3_300$Ivtot - list_siaya_3_300[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_siaya_3_300 <- malsim_odin[index,]

fit_siaya_3_300 <- my_sim_siaya_out_3_300_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_siaya_out_3_300, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_siaya_3_300, aes(x = t, y = Ivtot), col = "red")

#3x600
times_siaya_out_3_600_exp_decay <- my_sim_siaya_out_3_600_exp_decay %>%
  filter(between(t, endec_t_start, endec_t_end))

times_siaya_odin_3_600 <- my_sim_siaya_out_3_600 %>%
  filter(between(t, endec_t_start, endec_t_end))

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")

list_siaya_3_600 <- split(times_siaya_out_3_600_exp_decay, f = times_siaya_out_3_600_exp_decay$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_siaya_odin_3_600$Ivtot - list_siaya_3_600[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
best_fit_siaya_3_600 <- malsim_odin[index,]

fit_siaya_3_600 <- my_sim_siaya_out_3_600_exp_decay %>%
  filter(ref == index)

ggplot(my_sim_siaya_out_3_600, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = fit_siaya_3_600, aes(x = t, y = Ivtot), col = "red")

site <- "siaya"
dose <- c("1x400", "1x800", "3x300", "3x600")
endec_mu_vec <- c(best_fit_siaya_1_400$endec_mu, best_fit_siaya_1_800$endec_mu,
                  best_fit_siaya_3_300$endec_mu,
                  best_fit_siaya_3_600$endec_mu)

wane_vec <- c(best_fit_siaya_1_400$wane, best_fit_siaya_1_800$wane,
              best_fit_siaya_3_300$wane,
              best_fit_siaya_3_600$wane)

df_siaya <- data.frame(site, dose, endec_mu_vec, wane_vec)

df_list <- list(df_kwale, df_busia, df_hb, df_migori, df_siaya)

df_ivm_params <- do.call("rbind", df_list)
names(df_ivm_params) <- c("site", "dose", "mu_endec", "wane_endec")

df_kwale_params <- readRDS("analysis/IVER_ONE/kwale_param.rds") %>%
  mutate(site = "kwale")
df_busia_params <- readRDS("analysis/IVER_ONE/busia_param.rds") %>%
  mutate(site = "busia")
df_hb_params <- readRDS("analysis/IVER_ONE/homa_bay_param.rds") %>%
  mutate(site = "homa bay")
df_migori_params <- readRDS("analysis/IVER_ONE/migori_param.rds") %>%
  mutate(site = "migori")
df_siaya_params <- readRDS("analysis/IVER_ONE/siaya_param.rds") %>%
  mutate(site = "siaya")

list_sites <- list(df_kwale_params, df_busia_params, df_hb_params, df_migori_params, df_siaya_params)


df_params_all <- do.call("rbind", list_sites)

df_site_endec <- left_join(df_params_all, df_ivm_params, by = "site") %>%
  mutate(eff_len = case_when(
    dose == "1x400" ~ 58,
    TRUE ~ 60
  ))
saveRDS(df_site_endec, file = "analysis/IVER_ONE/output/df_site_endec.rds")

df_site_endec <-readRDS(file = "analysis/IVER_ONE/output/df_site_endec.rds")
