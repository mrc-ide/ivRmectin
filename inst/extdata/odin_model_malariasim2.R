


endec_killing <- mu + (endec_mu*avhc*ivm_cov)

deriv(Sv) <-  if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) -(FOIv*Sv) - endec_killing*Sv + betaa else - (FOIv*Sv) - mu*Sv + betaa
deriv(Ev[1]) <- if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) (FOIv*Sv) - Ev[1] - endec_killing*Ev[1] else (FOIv*Sv) - Ev[1] - mu**Ev[1]
deriv(Ev[2:10]) <-if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) Ev[i-1] - Ev[i] - endec_killing *Ev[i] else Ev[i-1] - Ev[i] - mu*Ev[i]
deriv(Iv) <- if (t >= (IVRM_sr)   && t < (IVRM_sr + eff_len)) Ev[10] - endec_killing*Iv else Ev[10] - mu*Iv


#do we need to add compartments for the ivm-killed mosquitoes and move between them because these could still transmit but just have a shorter lifespan?
#but drops in malariasim seemed okay, so may be okay
