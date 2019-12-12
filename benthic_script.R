library(tidyverse)
library(lubridate)
library(brms)
library(ggridges)
library(RCurl)
library(janitor)



# Load data ---------------------------------------------------------------
#benthic abundance (insects only)
ben_dm_tot <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/ben_dm_tot.csv")) 


#Bayesian brms model of benthic dry mass - download here or re-run using code below.
brm_ben_m2 <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/brm_ben_m2.RDS?raw=true"))





# Bayesian model ----------------------------------------------------------

#prior predictive - model based on priors only
# brm_ben_m2_prior <- brm(benthic_mg01~date*trt*taxon + (1|location), data=ben_dm_tot,family=Gamma(link="log"),
# prior=c(prior(normal(0,3),class="Intercept"),
# prior(normal(0,2),class="b"),
# prior(cauchy(0,1),class = "sd")), sample_prior = "only")
# 
# View prior predictive
# marginal_effects(brm_ben_m2_prior)

#full model
# brm_ben_m2 <- brm(benthic_mg01~date*trt*taxon + (1|location), data=ben_dm_tot,family=Gamma(link="log"),
# prior=c(prior(normal(0,3),class="Intercept"),
# prior(normal(0,2),class="b"),
# prior(cauchy(0,1),class = "sd")))

#check model and posterior predictive distributions (pp_check)
brm_ben_m2
pp_check(brm_ben_m2, type = "boxplot")

# saveRDS(brm_ben_m2, file = "brm_ben_m2.RDS")






# Extract conditional posteriors ------------------------------------------

date = unique(ben_dm_tot$date)
trt = unique(ben_dm_tot$trt)
taxon = unique(ben_dm_tot$taxon)
newdata = expand.grid(date,trt,taxon) %>% 
  rename(date = Var1,
         trt = Var2,
         taxon = Var3)

marg_post <- fitted(brm_ben_m2, summary = F, newdata = newdata, re_formula = NA)
names = newdata %>% unite("colnames", date:taxon)
colnames(marg_post) <- names$colnames

marg_post2 <- marg_post %>% as_tibble(marg_post)%>%
  mutate(iter = 1:nrow(marg_post)) %>% 
  gather(key,mg_dm_m2, -iter)%>%
  separate(key, c("date","trt","taxon"), sep = "_") %>% 
  mutate(date = mdy(as.factor(date)))

#total mg summed across all prey taxa
total_mg <- marg_post2 %>% 
  group_by(iter, date, trt) %>% 
  summarize(mg_dm = sum(mg_dm_m2)) %>% 
  mutate(group = "All insects")

#total mg chiro only
total_chiro <- marg_post2 %>%
  filter(taxon == "chironomidae") %>% 
  mutate(mg_dm = mg_dm_m2) %>% 
  mutate(group = "Chironomidae")











# Plot posterior ----------------------------------------------------------
# TOTAL MACROINVERT BIOMASS

# format raw data to plot on top of posterior
raw_ben_plot <- as_tibble(ben_dm_tot) %>% 
  mutate(date = mdy(as.factor(date)),
         trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, id, location, trt) %>% 
  summarize(mg_dm = sum(benthic_mg_dm_m2))



raw_chiro_plot <- as_tibble(ben_dm_tot) %>% 
  filter(taxon == "chironomidae") %>% 
  mutate(date = mdy(as.factor(date)),
         trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, id, location, trt) %>% 
  summarize(mg_dm = sum(benthic_mg_dm_m2))



#plot of raw data and posterior
plot_benthic <- total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                       "exc" = "no fish"))) %>%
  ggplot(aes(x = date, y = mg_dm, fill = trt)) +
  geom_boxplot(aes(group = interaction(trt, date)), outlier.shape = NA,
               position = position_dodge(width = 2),
               width = 1.5)+
  theme_classic()+
  scale_fill_grey(start = 0.9, end = 0.4)+
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  scale_y_continuous(limits = c(0, 4800),
                     breaks = c(0,1000, 2000, 3000, 4000, 5000)) +
  theme(legend.title = element_blank(),
        axis.title.x =element_blank(),
        text = element_text(size = 20))+
  geom_point(data = raw_ben_plot, aes(y = mg_dm, x = date, fill = trt), 
             position = position_dodge(width = 2),
             shape = 16, size = 1.3) +
  ylab(bquote('mg dry mass/'~m^2)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  #scale_y_log10()+
  ggtitle("b) Benthic invertebrates (84-91% chironomids)")

plot_benthic
ggsave(plot_benthic, file = "plot_benthic.tiff", dpi = 600, width = 7, height = 3.5, units = "in")
saveRDS(plot_benthic, file = "plot_benthic.rds")


# TOTAL CHIRONOMID BIOMASS
#plot of raw data dn posterior
plot_benthic_chiro <- total_chiro %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>%
  ggplot(aes(x = date, y = mg_dm_m2, fill = trt)) +
  geom_boxplot(aes(group = interaction(trt, date)), outlier.shape = NA,
               position = position_dodge(width = 2),
               width = 1.5)+
  theme_classic()+
  scale_fill_grey(start = 0.9, end = 0.4)+
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  scale_y_continuous(limits = c(0, 4800),
                     breaks = c(0,1000, 2000, 3000, 4000, 5000)) +
  theme(legend.title = element_blank(),
        axis.title.x =element_blank(),
        text = element_text(size = 20))+
  geom_point(data = raw_chiro_plot, aes(y = mg_dm, x = date, fill = trt), 
             position = position_dodge(width = 2),
             shape = 16, size = 1.3) +
  ylab(bquote('mg dry mass/'~m^2)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  #scale_y_log10()+
  ggtitle("b) Larval chironomids")

plot_benthic_chiro
ggsave(plot_benthic_chiro, file = "plot_benthic_chiro.tiff", dpi = 600, width = 7, height = 3.5, units = "in")
saveRDS(plot_benthic_chiro, file = "plot_benthic_chiro.rds")



# COMPARISON OF TOTAL EMERGENCE AND CHIRONOMID ONLY

plot_emerge_all_v_chiro <- total_chiro %>% 
  select(-taxon, -mg_dm_m2) %>% 
  bind_rows(total_mg) %>% 
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>%
  ggplot(aes(x = date, y = mg_dm, fill = trt)) +
  geom_boxplot(aes(group = interaction(trt, group, date)), outlier.shape = NA,
               position = position_dodge(width = 2),
               width = 1.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ group) +
  scale_fill_grey(start = 0.9, end = 0.4) +
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  scale_y_continuous(limits = c(0, 4800),
                     breaks = c(0,1000, 2000, 3000, 4000, 5000)) +
  theme(legend.title = element_blank(),
        axis.title.x =element_blank(),
        text = element_text(size = 10)) +
  # geom_point(data = raw_chiro_plot, aes(y = mg_dm, x = date, fill = trt), 
  #            position = position_dodge(width = 2),
  #            shape = 16, size = 1.3) +
  ylab(bquote('mg dry mass/'~m^2)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2) +
  #scale_y_log10()+
  NULL

plot_emerge_all_v_chiro
ggsave(plot_emerge_all_v_chiro, file = "plot_emerge_all_v_chiro.tiff", dpi = 600, width = 7, height = 3, units = "in")
saveRDS(plot_emerge_all_v_chiro, file = "plot_emerge_all_v_chiro.rds")




#posterior versus prior 
b_prior_ben = rnorm(4000, 0,2) #vector of priors from model

posts_prior_ben <- posterior_samples(brm_ben_m2) %>% 
  select(contains("b_")) %>% 
  gather(par, posterior) %>% 
  group_by(par) %>% 
  mutate(prior = ifelse(grepl("Intercept", par), rnorm(4000, 0,3), b_prior_ben)) %>% 
  gather(model, value, -par) %>% 
  ungroup() %>% 
  mutate(par = fct_relevel(as.factor(par), "b_Intercept", after = Inf))

plot_benthic_post_prior <- ggplot(posts_prior_ben, aes(x = value, y = par, fill = model)) +
  geom_density_ridges(quantile_lines = T, quantiles = 2) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(size = 8)) +
  geom_vline(xintercept = 0) +
  scale_fill_brewer(type = "qual", palette = 4) +
  ylab("Parameter") +
  xlab("Parameter value") +
  ggtitle("Benthic model prior versus posterior")

ggsave(plot_benthic_post_prior, file = "plot_benthic_post_prior.tiff", dpi = 600, width = 6.5, height = 7)
saveRDS(plot_benthic_post_prior, file = "plot_benthic_post_prior.rds")





# Summarize posterior -----------------------------------------------------

# TOTAL MACROINVERT BIOMASS
#summary stats by date and treatment
total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no_fish"))) %>% 
  group_by(date, trt) %>% 
  summarize(mean = mean(mg_dm),
            median = median(mg_dm),
            sd = sd(mg_dm),
            low95 = quantile(mg_dm, probs = 0.025),
            high95 = quantile(mg_dm, probs = 0.975))

#proportion difference by trt
total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, trt) %>% 
  pivot_wider(names_from = trt,
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(prop_reduction = fish/no_fish) %>% 
  summarize(mean = mean(prop_reduction),
            median = median(prop_reduction),
            sd = sd(prop_reduction),
            low95 = quantile(prop_reduction, probs = 0.025),
            high95 = quantile(prop_reduction, probs = 0.975))

#probability of a difference
total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, trt) %>% 
  pivot_wider(names_from = trt,
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(diff = fish-no_fish) %>% 
  summarize(prob_nofish_higher = sum(diff<0)/4000)


#proportion of benthics that are chironomids
total_chiro %>% mutate(chiro_mg = mg_dm) %>% 
  select(-mg_dm, -mg_dm_m2) %>% 
  left_join(total_mg) %>% 
  mutate(prop_chiro = chiro_mg/mg_dm) %>% 
  group_by(trt, date) %>% 
  summarize(mean = mean(prop_chiro),
            median = median(prop_chiro),
            sd = sd(prop_chiro),
            low95 = quantile(prop_chiro, probs = 0.025),
            high95 = quantile(prop_chiro, probs = 0.975))



---------------------------------------------------------------------------
# TOTAL CHIRONOMID BIOMASS

total_chiro %>%
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no_fish")),
         mg_dm = mg_dm_m2) %>% 
  select(trt, iter, mg_dm) %>% 
  group_by(trt) %>% 
  summarize(mean = mean(mg_dm),
            median = median(mg_dm),
            sd = sd(mg_dm),
            low95 = quantile(mg_dm, probs = 0.025),
            high95 = quantile(mg_dm, probs = 0.975))

#proportion difference by trt
total_chiro %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, trt) %>% 
  pivot_wider(names_from = trt,
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(prop_reduction = fish/no_fish) %>% 
  summarize(mean = mean(prop_reduction),
            median = median(prop_reduction),
            sd = sd(prop_reduction),
            low95 = quantile(prop_reduction, probs = 0.025),
            high95 = quantile(prop_reduction, probs = 0.975))

#probability of a difference
total_chiro %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, trt) %>% 
  select(-mg_dm_m2) %>% 
  pivot_wider(names_from = trt,
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(diff = fish-no_fish) %>% 
  summarize(prob_nofish_higher = sum(diff<0)/4000)





# Fish in cages -----------------------------------------------------------

fish_in_cages <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/fish_in_cages.csv")) %>% 
  mutate(date = mdy(date))

fish_in_cages_select <- fish_in_cages %>%  select(fish_in_cages, station)

benthic_plot_fishcage <- raw_ben_plot %>% 
  unite("station", c("id","location")) %>%
  left_join(fish_in_cages_select) %>% 
  mutate(fish_in_cages = ifelse(date == "2017-06-02",0, fish_in_cages)) %>% 
  ggplot(aes(x = date, y = mg_dm, fill = trt)) +
  geom_point(position = position_dodge(width = 2),
             shape = 21, alpha = 0.8, aes(size = fish_in_cages)) +
  scale_fill_grey(name = "")+
  scale_size_continuous(name = "# fish in cages on 2017-06-30")+
  theme_classic()+
  theme(axis.title.x =element_blank())+
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  ylab(bquote('mg dry mass/'~m^2)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  annotate("text",x=as.Date("2017-06-02")+7.5,y=3200,label="start of experiment")+
  geom_segment(aes(x = as.Date("2017-06-04"), y = 3200, 
                   xend=as.Date("2017-06-02"), yend = 3200),
               size = 0.1)+
  #scale_y_log10()+
  ggtitle("a) Benthic insects")+
  #annotate("text",x=as.Date("2017-06-02")+7.5,y=320,label="start of experiment")+
  #geom_segment(aes(x = as.Date("2017-06-05"), y = 320, xend=as.Date("2017-06-02"), yend = 320))+
  #scale_y_log10()+
  NULL



ggsave(benthic_plot_fishcage, file = "benthic_plot_fishcage.tiff", dpi = 600, width = 7, height = 3.5, units = "in")
saveRDS(benthic_plot_fishcage, file = "benthic_plot_fishcage.rds")

