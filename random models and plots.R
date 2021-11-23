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


site_data2 <- alldata %>% 
  group_by(site, PFT) %>% 
  summarize_all(mean) %>% 
  filter(site%in%sitenames) %>%
  ungroup() %>% 
  filter(PFT=="conifer")



sitenames <- unique(site_data$site)library(corrplot)
library(RColorBrewer)
M <-cor(dplyr::select(site_data2, c(LMA, VPDmean, VPDmax, precip, Tmean, Tmax)))
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))



plot_data <- alldata %>% 
  group_by(site, plot, PFT) %>% 
  summarize(mean.LMA=mean(LMA),
            n=n(),
            sd=sd(LMA),
            SE.LMA=sd(LMA)/sqrt(n),
            meanT=mean(Tmean),
            precip=mean(precip),
            lma.lower=mean.LMA -sd,
            lma.upper=mean.LMA +sd) %>% 
  filter(site%in%sitenames)

ggplot(site_data, aes(meanT, mean.LMA, color=PFT)) +
  geom_pointrange(aes(ymin=lma.lower, ymax=lma.upper)) +
  geom_smooth(method='lm')

hist(site_data$mean.LMA)
hist(log(site_data$mean.LMA))

library(lme4)

mod <- lmer(log(mean.LMA)~meanT + precip + (1|PFT), data=site_data)

library(brms)
library(tidybayes)

## spatial model----
bmod <- brm(mean.LMA|se(sd, sigma=TRUE) ~ 0 + Intercept + meanT*PFT, data=site_data, cores=4, iter=1000)
bmod2 <- brm(mean.LMA ~ 0 + Intercept + meanT*PFT, data=site_data, cores=4, iter=1000)
# bmod3 <- brm(mean.LMA ~ 0 + Intercept + meanT*PFT + (1|site), data=plot_data, cores=4, iter=1000)

beta <- bmod2 %>% 
  spread_draws(b_meanT, `b_meanT:PFTgrass`, `b_meanT:PFTconifer`) %>% 
  rename(b_broad=b_meanT, b_grass=`b_meanT:PFTgrass`, b_conif=`b_meanT:PFTconifer`) 

hist(beta$b_broad)


all_sub <- alldata %>% 
  filter(site%in%sitenames)

site_year <- alldata %>% 
  group_by(site) %>% 
  mutate(early=min(year)) %>% 
  group_by(site, PFT, year) %>% 
  summarize(mean.LMA=mean(LMA),
            meanT=mean(Tmean),
            precip=mean(precip),
            early=mean(early)) %>% 
  mutate(year2=year-(early-1)) %>% 
  filter(site%in%sitenames)

traittemp <- brm(mean.LMA ~ 0 + Intercept + meanT +(1|PFT + site), data=site_year, cores=4, iter=1000)
traittemp <- brm(mean.LMA ~ 0 + Intercept + meanT +(meanT|PFT/site), data=site_year, cores=4, iter=1000)

beta_time <- traittemp %>% 
  spread_draws(b_meanT, `b_meanT:PFTgrass`, `b_meanT:PFTconifer`) %>% 
  rename(b_broad=b_meanT, b_grass=`b_meanT:PFTgrass`, b_conif=`b_meanT:PFTconifer`) %>% 
  mutate(group="temporal")

beta_time <- traittemp %>% 
  spread_draws(b_meanT, r_PFT[pft, term], `r_PFT:site`[site, term]) %>% 
  mutate(b_site=b_meanT + r_PFT + `r_PFT:site`) %>%
  median_qi() %>% 
  ungroup() %>% 
  mutate(pft2=str_split_fixed(site, "_", 2)[,1]) %>% 
  filter(pft==pft2, term=="meanT")

beta_space <- bmod2 %>% 
  gather_draws(b_meanT, `b_meanT:PFTgrass`, `b_meanT:PFTconifer`) %>% 
  median_qi() %>% 
  mutate(pft=c("broadleaf", "conifer", "grass"))

p1 <- ggplot(beta_time, aes(y = site, x = b_site, xmin = b_site.lower, xmax = b_site.upper, color=pft)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0) +
  xlim(c(-2,1)) +
  facet_wrap(~pft, scales="free_y") +
  theme(legend.position = 'none') +
  ggtitle("temporal LMA slope")

p2 <- ggplot(beta_space, aes(y=pft, x = .value, xmin=.lower, xmax=.upper, color=pft)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0) +
  xlim(c(-2,1)) +
  facet_wrap(~pft, scales="free_y") +
  theme(legend.position = 'none') +
  ggtitle("spatial LMA slope")

allbetas <- left_join(beta_time, beta_space)

cowplot::plot_grid(p1, p2, ncol=1, rel_heights = c(1, 0.3), align = 'vh')

ggplot(allbetas, aes(y = site, x = b_site, xmin = b_site.lower, xmax = b_site.upper, color=pft)) +
  geom_pointinterval() +
  geom_vline(mapping = aes(xintercept=.value, color=pft)) +
  geom_vline(mapping = aes(xintercept=.lower, color=pft), linetype='dashed') +
  geom_vline(mapping = aes(xintercept=.upper, color=pft), linetype='dashed') +
  xlim(c(-2,1)) +
  facet_wrap(~pft, scales="free_y") +
  theme(legend.position = 'none') +
  ggtitle("temporal LMA slope") +
  theme_test() + 
  xlab("LMA slope estimate (LMA change per degree C)")
  



  

allbetas <- beta %>% 
  mutate(group="spatial") %>% 
  bind_rows(beta_time)

ggplot(allbetas, aes(b_broad, fill=group)) +
  geom_density(alpha=0.5) +
  ggtitle("Broadleaf LMA change per degree C")

ggplot(allbetas, aes(b_conif, fill=group)) +
  geom_density(alpha=0.5) +
  ggtitle("Conifer LMA change per degree C")

ggplot(allbetas, aes(b_grass, fill=group)) +
  geom_density(alpha=0.5) +
  ggtitle("Grassland LMA change per degree C")

ggplot(site_year, aes(year2, meanT, color=site)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  theme(legend.position = 'none')

starting <- site_year %>% 
  group_by(site) %>% 
  mutate(early=min(year)) %>% 
  filter(year==early)

yearchange <- site_year %>% 
  group_by(site) %>% 
  mutate(early=min(year),
         late=max(year),
         timediff=late-early) %>% 
  group_by(site, PFT) %>% 
  summarize(timediff=mean(timediff))

tempchange <- site_year %>% 
  group_by(site, PFT) %>% 
  mutate(early=min(year)+5,
         late=max(year)-5,
         timediff=late-early) %>% 
  mutate(group=ifelse(year<=early, "early", ifelse(year>=late, "late", NA))) %>% 
  filter(!is.na(group)) %>% 
  group_by(site, PFT, group, timediff) %>% 
  summarize(meanT=mean(meanT),) %>% 
  pivot_wider(id_cols=c("site", "PFT", "timediff"), names_from="group", values_from="meanT") %>% 
  mutate(tdiff=(late-early),
         tdiff_annual=tdiff/timediff)
  

bmod_time <- brm(mean.LMA ~ 0 + Intercept + meanT*PFT + (1|site), data=site_year, cores=4, iter=1000)


temp_change <- brm(meanT ~ 0 + Intercept + year + (year|site), data=site_year, cores=4, iter=1000)
temp_change2 <- brm(meanT ~ 0 + Intercept + year2, data=site_year, cores=4, iter=1000)

btemp <- fixef(temp_change2)[2,1]

tempchange <- tempchange %>% 
  mutate(tdiff_est=fixef(temp_change2)[2,1]*timediff)

beta2 <- beta %>% 
  pivot_longer(cols=c("b_broad", "b_grass", "b_conif"), names_to="PFT", values_to="beta")

beta2$PFT[beta2$PFT=="b_broad"] <- "broadleaf"
beta2$PFT[beta2$PFT=="b_grass"] <- "grass"
beta2$PFT[beta2$PFT=="b_conif"] <- "conifer"

simdata <- dplyr::select(starting, c(site, PFT, LMA_start=mean.LMA)) %>% 
  left_join(dplyr::select(tempchange, c(site, PFT, tdiff, tdiff_est, timediff))) %>% 
  left_join(beta2) %>% 
  mutate(LMA_end=LMA_start + beta*tdiff,
         est=beta)

sim_sum <- simdata %>% 
  group_by(site, PFT, LMA_start, timediff, tdiff_est, tdiff) %>% 
  summarize(med=median(est),
            lower=quantile(est, probs=0.025),
            upper=quantile(est, probs=0.975))

ggplot(sim_sum) +
  geom_pointrange(aes(x=med, y=site, xmin=lower, xmax=upper, color=PFT))+
  geom_vline(xintercept = 0, linetype='dashed')
  

lmachange <- site_year %>% 
  group_by(site, PFT) %>% 
  mutate(early=min(year)+5,
         late=max(year)-5,
         timediff=late-early) %>% 
  mutate(group=ifelse(year<=early, "early", ifelse(year>=late, "late", NA))) %>% 
  filter(!is.na(group)) %>% 
  group_by(site, PFT, group, timediff) %>% 
  summarize(LMA_mean=mean(mean.LMA),) %>% 
  pivot_wider(id_cols=c("site", "PFT", "timediff"), names_from="group", values_from="LMA_mean") %>% 
  mutate(lma_diff=(late-early),
         lma_diff_annual=lma_diff/timediff)



sim_with_dat <- sim_sum %>% 
  left_join(dplyr::select(lmachange, c(site, lma_diff))) %>% 
  mutate(lma_diff_per_deg=lma_diff/tdiff)

ggplot(sim_with_dat) +
  geom_pointrange(aes(x=med, y=site, xmin=lower, xmax=upper, color=PFT))+
  geom_point(aes(x=lma_diff, y=site)) +
  facet_wrap(~PFT, scales='free') +
  geom_vline(xintercept = 0, linetype='dashed')



site_year <- alldata %>% 
  group_by(site, PFT, year) %>% 
  summarize_if(.predicate = is.numeric, mean, na.rm=T) 


ggplot(plot_data, aes(Tmean, LMA, color=site)) +
  geom_point() +
  theme(legend.position = 'none')

ggplot(site_data, aes(Tmean, LMA, color=PFT)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(site_data, aes(precip, LMA, color=PFT)) +
  geom_point()+
  geom_smooth(method='lm')

ggplot(site_data, aes(VPDmean, LMA, color=PFT)) +
  geom_point()+
  geom_smooth(method='lm')

mod1 <- lm(LMA~ Tmean*PFT + precip*PFT, data=site_data)
summary(mod1)

mod_conf <- data.frame(confint(mod1))

diff_data <- diff_data %>% 
  mutate(delta.LMA.per.C=LMA.delta.mean/Tmean.delta)

ggplot(filter(diff_data, PFT=="conifer"), aes(delta.LMA.per.C)) +
  geom_histogram() + 
  geom_vline(xintercept = -2.903e-04, linetype='dashed') +
  geom_vline(xintercept = c(-6.829722e-04, 1.022746e-04), linetype='dashed', color='gray') 

ggplot(diff_data, aes(Tmean.delta, delta.LMA.per.C, color=PFT)) +
  geom_point() +
  geom_smooth(method = 'lm', se=F)

ggplot(diff_data, aes(precip.delta, fill=PFT)) + geom_density(alpha=0.4)
ggplot(diff_data, aes(VPDmean.delta, fill=PFT)) + geom_density(alpha=0.4)
ggplot(diff_data, aes(Tmean.delta, fill=PFT)) + geom_density(alpha=0.4)

  


