#load packages
source("src/packages.R")

# load data
d_tot <- read_csv("data/derived_data/d_tot.csv") # diet data
emerge_reu_mg <- read_csv("data/derived_data/emerge_reu_mg.csv") # emergence data
ben_dm_tot <- read_csv("data/derived_data/ben_dm_tot.csv") # benthic data
spider_abund <- read_csv("data/derived_data/spider_abund.csv") # spiderdata

# set options
rstan_options(autowrite = T)

# diet model --------------------------------------------------------------
# prior predictive

get_prior(g_dmtotal01 ~ 1 + date2*species*prey_stage, 
              family = Gamma(link = "log"),
              data = d_tot)

diet_brms_priors <- brm(g_dmtotal01 ~ 1 + date2*species*prey_stage, 
                        family = Gamma(link = "log"),
                        data = d_tot,
                        prior = c(prior(normal(-2,2), class = "Intercept"),
                                  prior(normal(0,2), class = "b")),
                 chains = 1, iter = 500, 
                 sample_prior = "only", #sample priors only
                 file = "models/prior_predictive/diet_brms_priors.rds",
                 file_refit = "on_change") 

prior_diet <- plot(conditional_effects(diet_brms_priors, effects = "date2:species", conditions = tibble(prey_stage = unique(d_tot$prey_stage))))
prior_diet$`date2:species` + scale_y_log10()
pp_check(diet_brms_priors, type = "boxplot")

# fit model
diet_brms <- brm(g_dmtotal01 ~ 1 + date2*species*prey_stage, 
                         family = Gamma(link = "log"),
                         data = d_tot,
                         prior = c(prior(normal(-2,2), class = "Intercept"),
                                   prior(normal(0,2), class = "b")),
                         cores = 4, 
                        file = "models/diet_brms.rds",
                        file_refit = "on_change")


#check model
pp_check(diet_brms, type = "boxplot") + scale_y_log10()




# emergence model ---------------------------------------------------------

#prior predictive
emerge_prior_dm_model <- brm(tot_mg_dm_m2_d ~ date*trt2 + (1|id), data=emerge_reu_mg,
                             family=Gamma(link="log"),
                             prior=c(prior(normal(1,2),class="Intercept"),
                                     prior(normal(0,4),class="b"),
                                     prior(exponential(1),class = "sd"),
                                     prior(gamma(2,1), class = "shape")),
                             sample_prior = "only",
                             chains = 1, iter = 500,
                             file = "models/prior_predictive/emerge_prior_dm_model.rds",
                             file_refit = "on_change")

emerge_reu_mg %>% 
  data_grid(date, trt2, id) %>% 
  add_predicted_draws(emerge_prior_dm_model) %>% 
  ggplot(aes(y = date, x= .prediction, color = trt2)) +
  geom_density_ridges() +
  scale_x_log10()


#full model
emerge_dm_model <- brm(tot_mg_dm_m2_d ~ date*trt2 + (1|id), data=emerge_reu_mg,
                             family=Gamma(link="log"),
                             prior=c(prior(normal(1,2),class="Intercept"),
                                     prior(normal(0,4),class="b"),
                                     prior(exponential(1),class = "sd"),
                                     prior(gamma(2,1), class = "shape")),
                             cores=4,
                       file = "models/emerge_dm_model.rds",
                       file_refit = "on_change")

plot(conditional_effects(emerge_dm_model, effects = "date:trt2"), points = T)
pp_check(emerge_dm_model, type = "boxplot") + scale_y_log10()

# benthic model -----------------------------------------------------------
#prior predictive

brm_ben_m2_prior <- brm(benthic_mg01~date*trt*taxon + (1|location),
                        data=ben_dm_tot,
                        family=Gamma(link="log"),
                        prior=c(prior(normal(5, 3),class="Intercept"),
                                prior(normal(0,1),class="b"),
                                prior(exponential(1), class = "sd"),
                                prior(gamma(1, 1), class = "shape")), 
                        chains = 1, iter = 500,
                        sample_prior = "only",
                        file = "models/prior_predictive/brm_ben_m2_prior.rds",
                        file_refit = "on_change")

conditional_effects(brm_ben_m2_prior, effects = "date:trt", conditions = tibble(taxon = unique(ben_dm_tot$taxon)))
pp_check(brm_ben_m2_prior, type = "boxplot") + scale_y_log10()

# full model
brm_ben_m2 <- update(brm_ben_m2_prior, 
                     cores = 4,
                     sample_prior = "no",
                     chains = 4, iter = 2000,
                     file = "models/brm_ben_m2.rds",
                     file_refit = "on_change")

#check model and posterior predictive distributions (pp_check)
test <- plot(conditional_effects(brm_ben_m2, effects = "date:trt", conditions = tibble(taxon = unique(ben_dm_tot$taxon))), points = T)
test$`date:trt` + scale_y_log10()
pp_check(brm_ben_m2, type = "boxplot") + scale_y_log10()


# spider model ------------------------------------------------------------

#emergence data for dates near spider collections
emerge_reu_mg_lasttwodates <- as_tibble(emerge_reu_mg) %>% 
  mutate(date2 = mdy(date)) %>% 
  filter(date2=="2017-06-16"|date2=="2017-06-23") %>% 
  mutate(spider_emerge_date = ifelse(date2=="2017-06-16",1,2),
         date3 = as.Date(date2),
         spider_emerge_date = as.Date(date3))%>% 
  select(date3,stat,stat2,loc,trt,tot_mg_dm_m2_d,spider_emerge_date) %>% 
  filter(trt!="Ambient")

#spider data to merge with emergence
spider_abund_tomerge <- as_tibble(spider_abund) %>% 
  mutate(spider_emerge_date = ymd(as.Date(date,"%m/%d/%Y")),
         date2 = ymd(as.Date(date,"%m/%d/%Y"))) %>% 
  select(date2,stat,stat2,loc,trt,spid,spider_emerge_date) %>% 
  filter(trt!="Ambient")

spider_emerge <- spider_abund_tomerge %>% 
  left_join(emerge_reu_mg_lasttwodates,by="stat2")



# prior predictive
spiders_prior_brm <- brm(spid~trt*date + (1|loc),data=spider_abund,
                         family=poisson(link="log"),
          prior=c(prior(normal(0,2),class="b"),
                  prior(normal(0,3),class="Intercept"),
                  prior(exponential(1), class = "sd")),
          iter = 1000, chains = 1, sample_prior = "only",
          file = "models/prior_predictive/spiders_prior_brm.rds",
          file_refit = "on_change")

conditional_effects(spiders_prior_brm, effects = "date:trt")
pp_check(spiders_prior_brm, type = "boxplot") + scale_y_log10()

spiders_brm <- brm(spid~trt*date + (1|loc),data=spider_abund,
                   family=poisson(link="log"),
      prior=c(prior(normal(0,2),class="b"),
             prior(normal(0,3),class="Intercept"),
            prior(exponential(1), class = "sd")),
   iter = 2000, chains = 4, cores=4,
   file = "models/spiders_brm.rds",
   file_refit = "on_change")

#check model
spiders_brm
pp_check(spiders_brm, type = "boxplot") + scale_y_log10()
