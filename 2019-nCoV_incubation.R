# 2019-nCOV_incubation.R
# The study approach is baesd on the paper 
# Backer et al. 2020. The incubation period of 2019-nCoV infections among travellers from Wuhan, China
# https://www.medrxiv.org/content/10.1101/2020.01.27.20018986v1

library(tidyverse)
library(rstan)
library(bayestestR)
library(loo)

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

# Download data from 
# https://docs.google.com/spreadsheets/d/1jS24DjSPVWa4iuxuD4OAXrE3QeI8c9BC1hSlqr-NMiU/edit#gid=1187587451
data <- read_tsv(file = "Line-list.tsv", skip = 1) %>%  
  dplyr::select(`reporting date`, location, gender, age, symptom_onset, 
         exposure_start, exposure_end, `visiting Wuhan`, `from Wuhan`) %>%
  filter(!is.na(exposure_end) & !is.na(symptom_onset))

data <- data %>% filter(`visiting Wuhan` == 1 | `from Wuhan` == 1)

tSymptomOnset <- as.integer((data$symptom_onset %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))
tStartExposure <- as.integer((data$exposure_start %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))
tEndExposure <- as.integer((data$exposure_end %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))

# exposure should ends before symptom onset
tEndExposure[tSymptomOnset < tEndExposure] <- tSymptomOnset[tSymptomOnset < tEndExposure]
# Constraint the minimal of tEndExposure and tSymptomOnset to 1
tSymptomOnset[which(tSymptomOnset - tEndExposure == 0)] <- tSymptomOnset[which(tSymptomOnset - tEndExposure == 0)]+1
# cases without maximum incubation period get start of exposure 21 days prior to first symptom onset
tStartExposure[is.na(tStartExposure)] <- min(tSymptomOnset) - 21

# input.data for stan analysis
input.data <- list(
  N = nrow(data),
  tStartExposure = tStartExposure,
  tEndExposure = tEndExposure,
  tSymptomOnset = tSymptomOnset)

# compile model
model <- stan(data = input.data, 
              chains = 0, 
              iter = 0,
              model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tSymptomOnset;
}

parameters{
  real<lower = 0> alphaInc; 	// Shape parameter of weibull distributed incubation period
  real<lower = 0> sigmaInc; 	// Scale parameter of weibull distributed incubation period
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
}

transformed parameters{
  vector[N] tE; 	// infection moment
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
}

model{
  // Contribution to likelihood of incubation period
  target += weibull_lpdf(tSymptomOnset -  tE  | alphaInc, sigmaInc);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = weibull_lpdf(tSymptomOnset[i] -  tE[i]  | alphaInc, sigmaInc);
  }
}
"
)

# Run simulation
stanfit <- stan(fit = model, data = input.data, 
                init = "random",
                iter = 10000, chains = 3)


# check results and convergence
print(stanfit)

# modelfit
LL = extract_log_lik(stanfit, parameter_name = 'log_lik')
loo(LL)

# Extract outputs
alpha <- rstan::extract(stanfit)$alphaInc %>% tail(500)
sigma <- rstan::extract(stanfit)$sigmaInc %>% tail(500)

#  Estimate 95% high density interval
hdi95 <- qweibull(0.95, shape=alpha, scale=sigma) %>% hdi(ci=0.9)

# Plot result
png(file = "fig1.png", width=2000, height=1500, res=300)
layout(matrix(c(1,2,2)), heights=c(1,2))
par(mar=c(0,4,2,2))
hist(qweibull(0.95, shape=alpha, scale=sigma), xlim=c(0,20), 
     main = "Estimated distribution of incubation days at 95th percentile",
     axes = F, col = "lightblue", ylab="")
y <- seq(0, 1, 0.01)
x <- qweibull(y, shape=alpha[1], scale=sigma[1])
par(mar=c(4,4,0,2))
plot(x, y, type="l", col="grey80", xlim=c(0,20), las=1,
     xlab="Incubation period (d)", ylab="CDF")
for(i in 2:500){
  x <- qweibull(y, shape=alpha[i], scale=sigma[i])
  lines(x, y, col="grey60")
}
abline(h=0.95, lty=2)
abline(v=14, lty=3, lwd=2, col="red")
lines(c(hdi95$CI_low, hdi95$CI_high), c(0.95,0.95), col = "lightblue", lwd=2)
text(14.1, 0, "Recommend quarantine days = 14", adj=0, col="red")
text(hdi95$CI_high, 0.9, "90% High density interval at 95th percentile", cex=1, adj=1, col="blue")
dev.off()
