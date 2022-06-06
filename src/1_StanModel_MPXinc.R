############################################################
#Purpose: Estimating the incubation period for monkeypox cases in the Netherlands, 2022
#Final edit: 6 June 2022
#Editor: Fumi Miura
############################################################
###Procedure
#0. install packages 
#1. data
#2. fit Stan models
#3. summary of model fits
#4. ggplots 
#5. references 
############################################################

###0. install packages -----
library(patchwork)
library(rstan)
library(loo)
library(tidyverse)#added
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
###1. data -----
raw_data <- read_csv("anonym_data.csv")

input_data <- list(N=length(raw_data$`Start date of exposure`), 
              tStartExposure=raw_data$`Start date of exposure`,
              tEndExposure=raw_data$`End date of exposure`,
              tSymptomOnset=raw_data$`Symptom onset`
)
###2. fit Stan models -----
##set distributions 
distributions <- c("weibull", "gamma", "lognormal")
code <- sprintf("
  data{
    int<lower=1> N;
    vector[N] tStartExposure;
    vector[N] tEndExposure;
    vector[N] tSymptomOnset;
  }
  parameters{
    real<lower=0> par[2];
    vector<lower=0, upper=1>[N] uE;	// Uniform value for sampling between start and end exposure
  }
  transformed parameters{
    vector[N] tE; 	// infection moment
    tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
  }
  model{
    // Contribution to likelihood of incubation period
    target += %s_lpdf(tSymptomOnset -  tE  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tSymptomOnset[i] -  tE[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
fit <- mapply(sampling, models, list(input_data), iter=10000, warmup=3000, chain=4)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

###3. summary of model fits -----
means <- cbind(pos$weibull[,2]*gamma(1+1/pos$weibull[,1]),
               pos$gamma[,1] / pos$gamma[,2],
               exp(pos$lognormal[,1]))
a_percentile <- c(0.025, 0.5, 0.975)
res <- apply(means, 2, quantile, a_percentile)
ll <- mapply(function(z) loo(extract_log_lik(z))$looic, fit)
waic <- mapply(function(z) waic(extract_log_lik(z))$waic, fit)
rbind(res, looIC=ll, WAIC=waic)

##Table-1 & Table-2
cens_w_percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qweibull(p = p, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025, 0.5, 0.975)))
colnames(cens_w_percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
cens_g_percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qgamma(p = p, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025, 0.5, 0.975)))
colnames(cens_g_percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
cens_ln_percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qlnorm(p = p, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025, 0.5, 0.975)))
colnames(cens_ln_percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
##Export Table-1 & Table-2
write.csv(round(cens_w_percentiles,1), "Weibull_percentile.csv")
write.csv(round(cens_g_percentiles,1), "Gamma_percentile.csv")
write.csv(round(cens_ln_percentiles,1), "Lognorm_percentile.csv")
write.csv(round(t(rbind(res, looIC=ll, WAIC=waic)),1), "meanICs.csv")

###4. ggplots -----
##make data frames for visualization
df <- data.frame(
  #Take mean values to draw emprical CDF
  inc_day = ((input_data$tSymptomOnset-input_data$tEndExposure)+(input_data$tSymptomOnset-input_data$tStartExposure))/2
)
x_plot <- seq(0,30,by=0.1)
Gam_plot <- as.data.frame(list(dose= x_plot, 
                               pred= sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.5))),
                               low = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025))),
                               upp = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.975)))
))
Wei_plot <- as.data.frame(list(dose= x_plot, 
                               pred= sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.5))),
                               low = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025))),
                               upp = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.975)))
))
ln_plot <- as.data.frame(list(dose= x_plot, 
                              pred= sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.5))),
                              low = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025))),
                              upp = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.975)))
))

##plot
gamma_ggplot <- ggplot(df, aes(x=inc_day)) +
  stat_ecdf(geom = "step")+ 
  xlim(c(0, 30))+
  geom_line(data=Gam_plot, aes(x=x_plot, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=Gam_plot, aes(x=x_plot,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 24)+
  labs(x="Incubation period (days)", y = "Proportion")+
  ggtitle("Gamma")

weibul_ggplot <- ggplot(df, aes(x=inc_day)) +
  stat_ecdf(geom = "step")+ 
  xlim(c(0, 30))+
  geom_line(data=Wei_plot, aes(x=x_plot, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=Wei_plot, aes(x=x_plot,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 24)+
  labs(x="Incubation period (days)", y = "Proportion")+
  ggtitle("Weibull")

lognorm_ggplot <- ggplot(df, aes(x=inc_day)) +
  stat_ecdf(geom = "step")+ 
  xlim(c(0, 30))+
  geom_line(data=ln_plot, aes(x=x_plot, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=ln_plot, aes(x=x_plot,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 24)+
  labs(x="Incubation period (days)", y = "Proportion")+
  ggtitle("Lognormal")

(lognorm_ggplot|gamma_ggplot|weibul_ggplot) + plot_annotation(tag_levels = 'A') #16inchi x 6 inchi

###5. references -----
#https://mc-stan.org/docs/2_29/functions-reference/lognormal.html
#http://mc-stan.org/rstanarm/reference/loo.stanreg.html
#https://mikuhatsune.hatenadiary.com/entry/2020/03/19/231812 (Japanese, to make three models at the same time)
