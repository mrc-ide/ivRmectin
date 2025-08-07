ivm_cov_par <- 0.7
ivm_min_age <-5
ivm_max_age <- 80
ivm_cov = ivm_cov_par*(exp(-ivm_min_age/21) - exp(-ivm_max_age/21))
prob_h_ivm <- 1/3

ivm_cov_in <- prob_h_ivm * ivm_cov
