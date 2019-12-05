# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(brms)
library(ggridges)
library(RCurl)
library(janitor)


# Load data ---------------------------------------------------------------

#full emergence abundance
emerge_reu_mg_long <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/emerge_reu_mg_long.csv"))
emerge_reu_mg <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/emerge_reu_mg.csv"))

#chiro_individual_lengths
chiro_ind <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/ind_ins.csv"))
#Bayesian model of chiro individual mass
chiro_mg_dm <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/emerge_chiro_mg_dm.RDS?raw=true"))
#Bayesian model of emergence (or re-run it below)
emerge_dm_model <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/emerge_dm_model.RDS?raw=true"))


# plot to check for cage-control vs ambient differences
emerge_reu_mg <- emerge_reu %>%
  replace(is.na(.), 0) %>%
  mutate(chiro_mg_dm = chir*.99,
         doli_mg_dm = doli*.65,
         tric_mg_dm = tric*.8,
         odo_mg_dm = odo*6.2,
         ephem_mg_dm = ephe*1.1,
         tot_mg_dm = chiro_mg_dm+doli_mg_dm+tric_mg_dm+odo_mg_dm+ephem_mg_dm,
         tot_mg_dm_m2_d = tot_mg_dm/area/days,
         date = mdy(date))

plot_emerge_cage_v_amb <- emerge_reu_mg %>% 
  filter(trt != "Exclusion",
         date != "2017-05-28") %>% 
  select(date, trt, chir,cera,doli,ephe,tric, simu, coleo, odo, hem) %>% 
  gather(taxon, abund, c(-date,-trt)) %>% 
  ggplot(aes(x = reorder(taxon, -abund), y = abund, color = trt, shape = trt))+
  geom_point(size = 2,position = position_jitterdodge(dodge.width = 0.6,
                                                      jitter.width = 0),
             alpha = .8)+
  theme_bw()+
  facet_grid(date~., scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 10))+
  scale_y_log10()+
  scale_color_grey(start = 0.2, end = 0.7)+
  xlab("prey_taxon")+
  ylab("Number per emergence trap")+
  #coord_flip() +
  NULL


ggsave(plot_emerge_cage_v_amb, file = "plot_emerge_cage_v_amb.tiff", dpi = 600, width = 6.5, height = 7)


# Length-mass regression - chiros -----------------------------------------

#plot to check for differences in chiro length over time and treatments
chiro_ind %>%
  mutate(date = mdy(date)) %>%
  ggplot(aes(x=date,y=mm,color=trt,group=interaction(trt,date)))+
  geom_violin(position=position_dodge(width=8))+
  geom_point(position=position_dodge(width=8),size=3,alpha=0.2)+
  #scale_y_log10()+
  theme_classic()+
  scale_color_grey()+
  ggtitle("Violin plot of chironomid lengths (mm)")+
  NULL

#add parameters (a,b) to convert length to mass
chiro_mg <- chiro_ind %>%
  mutate(a = 0.1, #parameter
         b = 1.57, #parameter
         param_source = "Sabo et al. 2002 (Table 2)", #source
         param_taxon = "Nematocera") %>% #source taxon
  mutate(mg_dm = a*mm^b) #generate biomass



# Fit model of chironomid average dry mass -----------------------------------------------------------

#chiro_mg_dm <- brm(mg_dm ~ 1, data=chiro_mg, family=Gamma(link="log"),
#                 prior=prior(normal(0,2),class="Intercept"))

#saveRDS(chiro_mg_dm, file = "emerge_chiro_mg_dm.RDS")

#check model and save
chiro_mg_dm
pp_check(chiro_mg_dm, type = "boxplot")
#saveRDS(chiro_mg_dm, file = "chiro_mg_dm.RDS")

#get mean and sd of chiro mass to sample from
post_chiro_mg <- posterior_samples(chiro_mg_dm)

chiro_adult_ind_shape <- mean(post_chiro_mg$shape)
chiro_adult_ind_scale <- mean(exp(post_chiro_mg$b_Intercept)/post_chiro_mg$shape)


#Convert emergence abundance to dry mass (mg)

#Chironomid dry mass based on model above using direct measures of mm
#Dry mass of other taxa were estimated from previous collections at the same site or nearby sites (Wesner et al. in review)

set.seed(2020)
emerge_reu_mg <- emerge_reu %>%
  replace(is.na(.), 0) %>%
  mutate(chiro_mg_dm = chir*.99,
         doli_mg_dm = doli*.65,
         tric_mg_dm = tric*.8,
         odo_mg_dm = odo*6.2,
         ephem_mg_dm = ephe*1.1,
         tot_mg_dm = chiro_mg_dm+doli_mg_dm+tric_mg_dm+odo_mg_dm+ephem_mg_dm,
         tot_mg_dm_m2_d = tot_mg_dm/area/days)
emerge_reu_mg$id <- paste(emerge_reu_mg$stat,emerge_reu_mg$loc)


emerge_reu_mg_long <- emerge_reu_mg %>% 
  select(date,id, trt2, days, area, 
         chiro_mg_dm, doli_mg_dm, 
         tric_mg_dm, odo_mg_dm, ephem_mg_dm) %>% 
  mutate(trt3 = ifelse(date == "05/28/17", "ctrl", trt2)) %>% 
  gather(taxon, mg_dm, c(-date, -id, -trt2, -days, -area))


# Bayesian model mg emergence ---------------------------------------------

#prior predictive
# emerge_prior_dm_model <- brm(tot_mg_dm_m2_d ~ date*trt2 + (1|id),data=emerge_reu_mg,family=Gamma(link="log"),
#                            prior=c(prior(normal(1,2),class="Intercept"),
#                                   prior=prior(normal(0,4),class="b")),
#                          sample_prior = "only")




#full model
# emerge_dm_model_taxon <- brm(mg_dm_m2_d01 ~ date*trt3*taxon + (1|id),data=emerge_reu_mg_long,
#                              family=Gamma(link="log"),
#                              prior=c(prior(normal(1,2),class="Intercept"),
#                                      prior(normal(0,4),class="b"),
#                                      prior(cauchy(0,1),class = "sd")),
#                              cores=4)

# saveRDS(emerge_dm_model_taxon, file = "emerge_dm_model.rds")

emerge_dm_model_taxon
pp_check(emerge_dm_model_taxon, type = "boxplot")

#make data to condition on
date <- unique(emerge_reu_mg_long$date)
trt3 <- unique(emerge_reu_mg_long$trt3)
taxon <- unique(emerge_reu_mg_long$taxon)
newdata <- expand.grid(date, trt2, taxon) %>% 
  rename(date = Var1,
         trt3 = Var2,
         taxon = Var3)

#extract posterior on each date an treatment
emerge_fit <- fitted(emerge_dm_model_taxon, newdata = newdata, summary=F, 
                     re_formula = NA)

names <- newdata %>% unite(colnames,c(date, trt3, taxon)) #names for the columns of the fitted estimates below
colnames(emerge_fit) <- names$colnames

post_emerge_mg_taxon <- as_tibble(emerge_fit) %>% 
  mutate(iter = 1:nrow(emerge_fit)) %>% 
  gather(key, mg_dm, -iter) %>% 
  separate(key, c("date","trt3","taxon"), sep = "_") 

post_emerge_mg_total <- post_emerge_mg_taxon %>% 
  group_by(date, trt3, iter) %>% 
  summarize(mg_dm = sum(mg_dm))

post_emerge_mg_chiro <- post_emerge_mg_taxon %>% 
  group_by(date, trt3, iter) %>% 
  filter(taxon == "chiro")

# Plot posteriors ---------------------------------------------------------
raw_data_plot <- emerge_reu_mg %>% 
  mutate(mg_dm = tot_mg_dm_m2_d,
         trt3 = if_else(date == "2017-05-28", "ctrl", as.character(trt2)),
         trt3 = str_replace_all(trt3,c("ctrl" = "fish",
                                      "exc" = "no fish")))

plot_emerge <- post_emerge_mg_total %>% ungroup() %>% 
  mutate(date = mdy(date),
         trt3 = str_replace_all(trt3,c("ctrl" = "fish",
                                  "exc" = "no fish"))) %>%
  ggplot(aes(x = date, y = mg_dm, fill = trt3)) +
  geom_boxplot(aes(group = interaction(trt3, date)), outlier.shape = NA,
               position = position_dodge(width = 2),
               width = 1.5) +
  theme_classic()+
  scale_fill_grey(start = 0.9, end = 0.4)+
  theme(legend.title = element_blank(),
        axis.title.x =element_blank(),
        text = element_text(size = 20))+
  coord_cartesian(ylim = c(0,2000),
                  xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  geom_point(data = raw_data_plot, aes(y = tot_mg_dm_m2_d,fill = trt3), 
             position = position_dodge(width = 2),
             shape = 16, size = 1.3) +
  annotate("text",x=as.Date("2017-06-02")+7.5,y=1500,label="start of experiment")+
  ylab(bquote('mg dry mass/'~m^2/d)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  geom_segment(aes(x = as.Date("2017-06-05"), y = 1500, xend=as.Date("2017-06-02"), yend = 1500))+
  #annotate("text",x=as.Date("2017-06-02")+7.5,y=320,label="start of experiment")+
  #geom_segment(aes(x = as.Date("2017-06-05"), y = 320, xend=as.Date("2017-06-02"), yend = 320))+
  #scale_y_log10()+
  ggtitle("a) Emerging insects (>99% chironomids)")


ggsave(plot_emerge, file = "plot_emerge.tiff", dpi = 600, width = 7, height = 3.5, units = "in")
saveRDS(plot_emerge, file = "plot_emerge.rds")


# Summarize posteriors ----------------------------------------------------

#summary of each date and treatment
post_emerge_mg %>% 
  mutate(date = mdy(date),
         trt2 = str_replace_all(trt2,c("ctrl" = "no fish",
                                       "exc" = "fish"))) %>% 
  group_by(date, trt2) %>% 
  summarize(mean = mean(mg_dm),
            median = median(mg_dm),
            sd = sd(mg_dm),
            low95 = quantile(mg_dm, probs = 0.025),
            high95 = quantile(mg_dm, probs = 0.975)) 


#summary of cumulative emergence after experiment began
sum_emerge <- post_emerge_mg %>% 
  mutate(date = mdy(date),
         trt2 = str_replace_all(trt2,c("ctrl" = "no fish",
                                       "exc" = "fish"))) %>% 
  pivot_wider(names_from = date, 
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(tot_emerge = x2017_06_06 + x2017_06_12 + x2017_06_16 + x2017_06_23) %>% 
  select(iter, trt2, tot_emerge) %>% 
  group_by(trt2) %>% 
  summarize(mean = mean(tot_emerge),
            median = median(tot_emerge),
            sd = sd(tot_emerge),
            low95 = quantile(tot_emerge, probs = 0.025),
            high95 = quantile(tot_emerge, probs = 0.975))

write.csv(sum_emerge, file = "sum_emerge.csv")  



#cumulative difference
post_emerge_mg %>% 
  mutate(date = mdy(date),
         trt2 = str_replace_all(trt2,c("ctrl" = "no fish",
                                       "exc" = "fish"))) %>% 
  pivot_wider(names_from = date, 
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(tot_emerge = x2017_06_06 + x2017_06_12 + x2017_06_16 + x2017_06_23) %>% 
  select(iter, trt2, tot_emerge) %>% 
  group_by(trt2) %>% 
  pivot_wider(names_from = trt2, values_from = tot_emerge) %>% 
  clean_names() %>% 
  mutate(prop_reduce = 1 - no_fish/fish,
         diff = fish - no_fish) %>% 
  summarize(mean = mean(prop_reduce),
            median = median(prop_reduce),
            sd = sd(prop_reduce),
            low95 = quantile(prop_reduce, probs = 0.025),
            high95 = quantile(prop_reduce, probs = 0.975),
            mean_diff = mean(diff),
            median_diff = median(diff),
            sd_diff = sd(diff),
            low95_diff = quantile(diff, probs = 0.025),
            high95_diff = quantile(diff, probs = 0.975))

#how much more emergence over 28 days
post_emerge_mg %>% 
  mutate(date = mdy(date),
         trt2 = str_replace_all(trt2,c("ctrl" = "no fish",
                                       "exc" = "fish"))) %>% 
  pivot_wider(names_from = date, 
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(tot_emerge = x2017_06_06 + x2017_06_12 + x2017_06_16 + x2017_06_23) %>% 
  select(iter, trt2, tot_emerge) %>% 
  group_by(trt2) %>%
  pivot_wider(names_from = trt2, values_from = tot_emerge) %>% 
  clean_names() %>% 
  mutate(daily_diff = fish - no_fish,
         diff_28d = daily_diff*28) %>% 
  summarize(mean = mean(diff_28d),
            median = median(diff_28d),
            sd = sd(diff_28d),
            low95 = quantile(diff_28d, probs = 0.025),
            high95 = quantile(diff_28d, probs = 0.975))

#proportion_chiros
post_emerge_mg_chiro %>% 
  mutate(mg_dm_chiro = mg_dm) %>% 
  select(-mg_dm, -taxon) %>% 
  left_join(post_emerge_mg_total) %>% 
  mutate(prop_chiro = mg_dm_chiro/mg_dm) %>% 
  group_by(date,trt3) %>% 
  summarize(mean = mean(prop_chiro),
            median = median(prop_chiro),
            sd = sd(prop_chiro),
            low95 = quantile(prop_chiro, probs = 0.025),
            high95 = quantile(prop_chiro, probs = 0.975))




# Fish in cages -----------------------------------------------------------

#Fish were found in 5 of the exclusion cages at the end of the experiment. To check their
#influence, plot emergence with fish in cages as an indicator. Did having fish tend to reduce
#emergence in the cages?

#get data on fish in cages at the end of the experiment
fish_in_cages <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/fish_in_cages.csv")) %>% 
  mutate(date = mdy(date))

fish_in_cages_select <- fish_in_cages %>%  select(fish_in_cages, station)


#plot it
emerge_plot_fishcage <- raw_data_plot %>% 
  unite("station", c("stat","loc")) %>%
  left_join(fish_in_cages_select) %>% 
  mutate(fish_in_cages = ifelse(date == "2017-05-28",0,fish_in_cages),
         fish_01 = ifelse(fish_in_cages == 0,0, 1)) %>% 
  ggplot(aes(x = date, y = mg_dm, fill = trt2, size = fish_in_cages)) +
  geom_point(aes(y = tot_mg_dm_m2_d), 
             position = position_dodge(width = 2),
             shape = 21, alpha = 0.8) +
  #geom_text_repel(aes(label = station))+
  scale_fill_grey(name = "")+
  scale_size_continuous(name = "# fish in cages on 2017-06-30")+
  theme_classic()+
  theme(axis.title.x =element_blank())+
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  ylab(bquote('mg dry mass/'~m^2/d)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  #annotate("text",x=as.Date("2017-06-02")+7.5,y=320,label="start of experiment")+
  #geom_segment(aes(x = as.Date("2017-06-05"), y = 320, xend=as.Date("2017-06-02"), yend = 320))+
  #scale_y_log10()+
  ggtitle("b) Emerging aquatic insects") +
  NULL

ggsave(emerge_plot_fishcage, file = "emerge_plot_fishcage.tiff", dpi = 600, width = 7, height = 3.5, units = "in")
saveRDS(emerge_plot_fishcage, file = "emerge_plot_fishcage.rds")

