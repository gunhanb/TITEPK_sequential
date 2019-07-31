#####################################################################
### 31.07.2019
### Authors: Burak Kürsad Günhan and Sebastian Weber
#####################################################################
### Code to implement TITE-PK method
### This is the accompanying code to paper 
### "Sequential phase I dose-escalation trials with multiple schedules"
### See the paper for details
#####################################################################


##############################################################################
### LICENCE
### 
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
### AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
### OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
### SOFTWARE.
################################################################################




########################
## Load the R packages
#######################
library(rstan)         ## Bayesian estimation
library(bayesplot)     
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(rcpp.cache.dir="cache")
library(ggplot2)
library(gridExtra)
library(functional)
library(plyr)
library(assertthat)
library(loo)          ## Model comparison

########################
### Settings for Plots
#########################
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
             strip.text.x = element_text(size = 12),
             axis.text = element_text(size = 9), 
             axis.title = element_text(size = 9, face = "bold"))

############################################
## Load the R files for TITE-PK method
source("src/utils.R")
source("src/utils_tite_pk.R")
source("lib/tools.R")
###########################################

############################################
## Compiling Stan files for TITE-PK method
###########################################
stan_model_code_titepk <- stanc_builder("src/tite_pk.stan")
cat(stan_model_code_titepk$model_code, file="src/tite_pk_run.stan")
cat("\n", file="src/tite_pk_run.stan", append = TRUE)
## compile Stan model functions via stanc_builder which avoids model
## name obfuscation and then admits caching in case nothing changes
model_exposed_cached <- stanc_builder("src/tite_pk_run.stan")
cat("Exposing Stan functions...")
expose_stan_functions(model_exposed_cached)
cat("done.\n")
## Stan code for estimation
model_titepk <- stan_model("src/tite_pk_run.stan", "TITE-PK One parameter model")

############################################
## Fixed Pseudo-PK model parameters
## T_e: half elimination rate constant
## Calculated from PK analysis (mean of the estimated T_e values)
T_e = log(2) / 30
## k_eff: kinteic constant which govern delay betwenn conc. in central compart and effect site
## calculated from PK analysis
k_e = exp(0.41)
#########################################################
## Setting up the TITE-PK Regimens to be considered
#########################################################
## This is similar to a NONMEM data structure
## Choosing the reference regimen
ref_dose    = 5
ref_tau     = 24
## Since cycle length is given as 21 days
tref_month  = 0.75
## Firstly consider the Daily RAD dataset
tau_drug    = 24
addl        = 730
# calculating referecne scale
source("src/pk_param.R")
# Creating provisional regimen (Doses and frequencies)
source("src/regimen_ref.R")
regimen_ref_packed
dosings <- regimen_ref_packed
# provisional Doses in the trial
Doses_daily   <- c(2.5, 5, 7.5, 10)
Doses   <- Doses_daily 
dosings <- dosings[rep(seq_len(nrow(dosings)), each=length(Doses)),]
dosings$amt  <- Doses
dosings$lamt <- log(Doses)
regimens_daily = dosings
pk_calc <- Curry(pk_model, theta=params$theta, lscale=ref_lscale)
# calculating C(t_ref), E(t_ref), AUC_E(t_ref) for reference regimen
pksys <- ddply(dosings, .(amt), pk_subject, time = 1:tref_h, pk_fun = pk_calc)

## TITE-PK: Prior value for log(beta) paramater
Prior_titepk = c(cloglog(0.175), 1.25)

###########################################################
## TITE-PK
############################################################
## Stan model for estimatiosn
## Data needs to be in a list format for Stan
## Load the dataset
# For our model
dat_RAD_daily <- read.csv("data/dat_RAD_daily")
dat_RAD_daily$X <- NULL
head(dat_RAD_daily)
time_unit <- 24 * 7 * 4
stanDat_posteriors_daily <- c(dat_RAD_daily,
                              list(N = nrow(dat_RAD_daily), 
                                   theta = log(c(T_e, k_e)),
                                   ref_dose = ref_dose,
                                   ref_tau = ref_tau,
                                   tref_month = tref_month,
                                   Nregimens = nrow(regimens_daily),
                                   regimens_lamt = regimens_daily$lamt,
                                   regimens_tau = regimens_daily$tau,
                                   addls = regimens_daily[, 7],
                                   params_prior =  Prior_titepk))



### Calculate aprior DLT probs
system.time(fit_posteriors_daily_tte <- sampling(model_titepk
                                                 ,stanDat_posteriors_daily
                                                 ,cores=1
                                                 ,init=1
                                                 ,chains=4
                                                 ,seed=12313
                                                 ,refresh=0
                                                 ,iter=2000
                                                 ,control=list(adapt_delta=0.975)
                                                 ,open_progress=FALSE))
## Timing
#user  system elapsed 
#0.783   0.029   0.908 
## Convergence diagnostics
#library(shinystan)
#launch_shinystan(fit_posteriors_daily_tte)
## Extract the necessary values
prefs_posteriors_daily_tte <- matrix(summary(fit_posteriors_daily_tte)$summary[c("P_dose[1]", "P_dose[2]",
                                                                                 "P_dose[3]", "P_dose[4]"), 
                                                                               c(4, 5, 6, 7, 8)], ncol = 5)


#######################################
## Plot the posterior: daily Figure 2B
d <- data.frame(x = factor(c(2.5, 5, 7.5, 10)),
                y = as.vector(rbind(prefs_posteriors_daily_tte[, 3])),
                ylo = as.vector(rbind(prefs_posteriors_daily_tte[, 1])),
                yhi =  as.vector(rbind(prefs_posteriors_daily_tte[, 5])),
                yinlo = as.vector(rbind(prefs_posteriors_daily_tte[, 2])),
                yinhi = as.vector(rbind(prefs_posteriors_daily_tte[, 4])))

args_inner <- list(mapping = aes_(ymax = ~ yinhi, ymin = ~ yinlo, x = ~ x),
                   size = 1.5,
                   show.legend = FALSE,
                   position = position_dodge(width = 0.50))
layer_inner <- do.call(geom_linerange, args_inner)


plot.dlts.posteriors.daily <- ggplot(d, aes(x = x, y = y, ymin = ylo, ymax = yhi)) + 
  geom_pointrange(position = position_dodge(width = 0.50),
                  aes()) + 
  layer_inner +
  geom_hline(aes(yintercept=0.16), lty = 2) +
  geom_hline(aes(yintercept=0.33), lty = 2) +
  xlab("Doses (mg/daily)") +
  ylab("DLT probs.") +
  ggtitle("B (Posterior: Daily)") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
  theme(legend.title.align=0.5) 


plot.dlts.posteriors.daily

#############
# PRIORS
#############
### One parameter model
## Stan model for prior derivations
stan_model_code_titepk_priors <- stanc_builder("src/tite_pk_priors.stan")
cat(stan_model_code_titepk_priors$model_code, file="src/tite_pk_priors_run.stan")
cat("\n", file="src/tite_pk_priors_run.stan", append = TRUE)
model_priors_titepk <- stan_model("src/tite_pk_priors_run.stan", "TITE-PK One parameter PriorsOnly model")
## Data needs to be in a list format for Stan
stanDat_priors_daily <-  list(theta = log(c(T_e, k_e)),
                              ref_dose = ref_dose,
                              ref_tau = ref_tau,
                              tref_month = tref_month,
                              regimens_lamt = regimens_daily$lamt,
                              regimens_tau = regimens_daily$tau,
                              Nregimens = nrow(regimens_daily),
                              addls = regimens_daily[, 7],
                              params_prior =  Prior_titepk)

### Calculate aprior DLT probs
fit_priors_daily_tte <- sampling(model_priors_titepk
                                 ,stanDat_priors_daily
                                 ,cores=1
                                 ,init=1
                                 ,chains=1
                                 ,seed=12313
                                 ,refresh=0
                                 ,control=list(adapt_delta=0.975)
                                 ,open_progress=FALSE,
                                 algorithm = "Fixed_param")
## Extract the necessary values
prefs_priors_daily_tte <- matrix(summary(fit_priors_daily_tte)$summary[c("P_dose[1]", "P_dose[2]",
                                                                         "P_dose[3]", "P_dose[4]"), 
                                                                       c(4, 5, 6, 7, 8)], ncol = 5)



#################################
### Figure 2A
#################################
d <- data.frame(x = factor(c(2.5, 5, 7.5, 10)),
                y = as.vector(rbind(prefs_priors_daily_tte[, 3])),
                ylo = as.vector(rbind(prefs_priors_daily_tte[, 1])),
                yhi =  as.vector(rbind(prefs_priors_daily_tte[, 5])),
                yinlo = as.vector(rbind(prefs_priors_daily_tte[, 2])),
                yinhi = as.vector(rbind(prefs_priors_daily_tte[, 4])))


args_inner <- list(mapping = aes_(ymax = ~ yinhi, ymin = ~ yinlo, x = ~ x),
                   size = 1.5,
                   show.legend = FALSE,
                   position = position_dodge(width = 0.50))
layer_inner <- do.call(geom_linerange, args_inner)


plot.dlts.priors.daily <- ggplot(d, aes(x = x, y = y, ymin = ylo, ymax = yhi)) + 
  geom_pointrange(position = position_dodge(width = 0.50),
                  aes()) + 
  layer_inner +
  geom_hline(aes(yintercept=0.16), lty = 2) +
  geom_hline(aes(yintercept=0.33), lty = 2) +
  xlab("Doses (mg/daily)") +
  ylab("DLT probs.") +
  ggtitle("A (Prior)") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
  theme(legend.title.align=0.5) 

plot.dlts.priors.daily


######
###### Considering both daily and weekly dosing regimsens
######
############################################
## Fixed Pseudo-PK model parameters
## T_e: half elimination rate constant
## fixed PK parameters

## WEEKLY DOSING
tau = 24 * 7
addl = 120
# Calculate refernce scale
source("src/pk_param.R")
pk_calc <- Curry(pk_model, theta=params$theta, lscale=ref_lscale)
tau_drug <- tau
# Calculating aucs
source("src/regimen_ref.R")
Doses_weekly <- c(20, 30, 50)
dosings_weekly <- regimen_ref_packed
dosings_weekly <- dosings_weekly[rep(seq_len(nrow(dosings_weekly)), each=length(Doses_weekly)),]
dosings_weekly$amt  <- Doses_weekly
dosings_weekly$lamt <- log(Doses_weekly)
pksys <- ddply(dosings_weekly, .(amt), pk_subject, time = tref_h, pk_fun = pk_calc)
# Only taking log(AUC_E(tref)) values for calculating dose specific DLT rates
laucs_weekly <- pksys[, c(4)]

regimens <- rbind(dosings, dosings_weekly)
## Note that "dosings" is only for daily dosing

## combined analysis
## Load the dataset
# For our model
dat_RAD_combined <- read.csv("data/dat_RAD_combined")
dat_RAD_combined$X <- NULL
head(dat_RAD_combined)
time_unit <- 24 * 7 * 4
stanDat_combined <- c(dat_RAD_combined,
                      list(N = nrow(dat_RAD_combined), 
                           theta = log(c(T_e, k_e)),
                           ref_dose = ref_dose,
                           ref_tau = ref_tau,
                           tref_month = tref_month,
                           regimens_lamt = regimens$lamt,
                           regimens_tau = regimens$tau,
                           Nregimens = nrow(regimens),
                           addls = regimens[, 7],
                           params_prior =  Prior_titepk))

###
# posterior
##
fit_posterior_combined_tte <- sampling(model_titepk
                                       ,stanDat_combined
                                       ,cores=1
                                       ,init=1
                                       ,chains=4
                                       ,seed=12313
                                       ,refresh=0
                                       ,iter=2000
                                       ,control=list(adapt_delta=0.975)
                                       ,open_progress=FALSE)
prefs_posterior_combined_tte <- matrix(summary(fit_posterior_combined_tte)$summary[c("P_dose[1]", "P_dose[2]",
                                                                                     "P_dose[3]", "P_dose[4]"), 4:8], ncol = 5)




###
### Figure 2C
###
d <- data.frame(x = factor(c(2.5, 5, 7.5, 10)),
                y = as.vector(rbind(prefs_posterior_combined_tte[, 3])),
                ylo = as.vector(rbind(prefs_posterior_combined_tte[, 1])),
                yhi =  as.vector(rbind(prefs_posterior_combined_tte[, 5])),
                yinlo = as.vector(rbind(prefs_posterior_combined_tte[, 2])),
                yinhi = as.vector(rbind(prefs_posterior_combined_tte[, 4])))


args_inner <- list(mapping = aes_(ymax = ~ yinhi, ymin = ~ yinlo, x = ~ x),
                   size = 1.5,
                   show.legend = FALSE,
                   position = position_dodge(width = 0.50))

layer_inner <- do.call(geom_linerange, args_inner)


plot.dlts.posteriors.combined <- ggplot(d, aes(x = x, y = y, ymin = ylo, ymax = yhi)) + 
  geom_pointrange(position = position_dodge(width = 0.50),
                  aes()) + 
  layer_inner +
  geom_hline(aes(yintercept=0.16), lty = 2) +
  geom_hline(aes(yintercept=0.33), lty = 2) +
  xlab("Doses (mg/daily)") +
  ylab("DLT probs.") +
  ggtitle("C (Posterior: Daily + Weekly)") + 
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
  theme(legend.title.align=0.5) 


plot.dlts.posteriors.combined
## Plot all figures together
grid.arrange(plot.dlts.priors.daily, plot.dlts.posteriors.daily, 
             plot.dlts.posteriors.combined, ncol = 3)

############
# END
############
