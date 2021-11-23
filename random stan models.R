library(tidyverse)

alldata <- read.csv("LMA_and_climate_allData.csv") %>% 
  mutate(plot=ifelse(!is.na(plotID), plotID, plot),
         LMA=LMA*1000) %>% 
  group_by(site, plot) %>% 
  mutate(nobs=n()) %>% 
  filter_at(vars(year, site, LMA, PFT), all_vars(!is.na(.)))
# filter(nobs>=10)


site_data <- alldata %>% 
  group_by(site, PFT) %>% 
  summarize(mean.LMA=mean(LMA),
            n=n(),
            sd=sd(LMA),
            SE.LMA=sd(LMA)/sqrt(n),
            meanT=mean(Tmean),
            precip=mean(precip),
            lma.lower=mean.LMA -sd,
            lma.upper=mean.LMA +sd) %>% 
  filter(n>5) %>% 
  filter(PFT!="mixed")

sitenames <- unique(site_data$site)

plot_data <- alldata %>% 
  group_by(site) %>% 
  mutate(plot2=as.integer(as.factor(plot)),
         plot2=ifelse(is.na(plot2), 1, plot2)) %>% 
  group_by(site, plot, plot2, PFT) %>% 
  summarize(mean.LMA=mean(LMA),
            n=n(),
            sd=sd(LMA),
            SE.LMA=sd(LMA)/sqrt(n),
            meanT=mean(Tmean),
            precip=mean(precip),
            lma.lower=mean.LMA -sd,
            lma.upper=mean.LMA +sd) %>% 
  group_by(site, PFT) %>% 
  mutate(meanT=mean(meanT)) %>% 
  filter(site%in%sitenames)

spatial_temp <- "
data {
  int<lower=0> J; 
  int<lower=0> N; 
  int<lower=1,upper=J> site[N];
  vector[N] temp;
  vector[N] y;
} 
parameters {
  vector[J] a;
  vector[1] b;
  real mu_a;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_y;
} 
transformed parameters {
  vector[N] y_hat;
  vector[N] m;

  for (i in 1:N) {
    m[i] = a[site[i]] + temp[i] * b[1];
    y_hat[i] = m[i];
  }
}
model {
  mu_a ~ normal(0, 1);
  a ~ normal(mu_a, sigma_a);
  b ~ normal(0, 1);
  y ~ normal(y_hat, sigma_y);
}

generated quantities {
real y_rep[N] = normal_rng(y_hat, sigma_y);
}"

data2 <- list(J=length(sitenames), N=nrow(plot_data), site=as.integer(factor(plot_data$site)),
              y=log(plot_data$mean.LMA),
              temp=plot_data$meanT)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


fit<- stan(model_code=spatial_temp, data=data2, iter=1000, warmup=500, chains=4, cores=4)
posterior <- rstan::extract(fit, inc_warmup = TRUE, permuted = FALSE)

plot(fit)

y_rep <- rstan::extract(fit, inc_warmup = TRUE, permuted = FALSE, pars="y_rep")

library(tidybayes)
library(shinystan)

ppc_dens_overlay(data2$y, y_rep[1:50, 2,])

bfit <- brm(log(mean.LMA) ~ meanT + (1|site) + (meanT|PFT), data=plot_data, prior=prior1, iter=1000, cores=4)

bfit %>% 
  spread_draws(b_Intercept, r_site[site], r_PFT[PFT, term]) %>%
  median_qi(site_mean=b_Intercept + r_site + r_PFT) %>% 
  filter(term=="Intercept") %>% 
  left_join(dplyr::select(site_data, c(site, meanT, PFT))) %>% 
  filter(!is.na(meanT)) %>% 
  ggplot(aes(y = reorder(site, meanT), x = site_mean, xmin=.lower, xmax=.upper, color=PFT)) +
    geom_pointinterval()

bfit %>% 
  spread_draws(b_meanT, r_PFT[PFT, term]) %>%
  median_qi(pft_slope=b_meanT + r_PFT) %>% 
  filter(term=="meanT") %>% 
  ggplot(aes(y = PFT, x = pft_slope, xmin=.lower, xmax=.upper, color=PFT)) +
  geom_pointinterval()


bfit %>% 
  spread_draws(b_Intercept, r_site[site]) %>%
  median_qi(site_mean=b_Intercept + r_site) %>% 
  left_join(dplyr::select(site_data, c(site, meanT, PFT))) %>% 
  ggplot(aes(y = reorder(site, meanT), x = site_mean, xmin=.lower, xmax=.upper, color=PFT)) +
  geom_pointinterval()

