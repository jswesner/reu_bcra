library(tidyverse)
library(lubridate)
library(brms)
library(ggridges)
library(RCurl)
library(janitor)



# Load data ---------------------------------------------------------------
#benthic abundance (insects only)
ben <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/benthic62617.csv")) %>% 
  clean_names() %>% 
  mutate(chironomidae = chiro_larvae + pupae) %>% 
  select(-chiro_larvae, -pupae)

#benthic dry mass of individuals
ben_dm_ind <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/benthic_drymass_141517.csv")) %>% 
  filter(mg_dm_per_ind > 0) %>% 
  mutate(taxon = str_replace(taxon, "Baetidae", "ephemeroptera"),
         taxon = str_replace(taxon, "Caenidae", "ephemeroptera"))

#Bayesian brms model of length-weight regression
brm_ind_mg <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/benthic_ind_mg.RDS?raw=true"))

#Bayesian model of benthic individual dry mass
#brm_ind_mg <- brm(mg_dm_per_ind ~ 0 + taxon + (1|experiment), data = ben_dm_ind,
#                  family = Gamma(link = "log"),
#                  prior = c(prior(normal(0,1), class = "b")),
#                 cores = 4)

#check model
# brm_ind_mg
# pp_check(brm_ind_mg, type = "boxplot")

#saveRDS(brm_ind_mg, file = "benthic_ind_mg.RDS")

#extract posterior and estimate mean
ind_mg_post <- posterior_samples(brm_ind_mg) %>% 
  select(starts_with('b')) %>% 
  gather()

ind_mg_mean <- ind_mg_post %>% 
  group_by(key) %>% 
  summarize(mean_mgdm = mean(exp(value))) %>% 
  mutate(taxon = tolower(str_remove(key, "b_taxon"))) %>% 
  select(-key) %>% 
  add_row(taxon="tipulidae",mean_mgdm=1.4) %>% 
  mutate(taxon = ifelse(taxon == "worms", "olgigochaetes", taxon))
#replace Caenidae with Ephem. #most mayflies are Caenids anyway at this 
#site, so this allows easer matching of names later on.
#Tipulidae were not measured for length. Instead, the mean of tipulid 
#lengths and regression to mass from Benke et al. (1999) was used to 
#estimate Tipulid length. Thus, there is no variation around this estimate
#(or at least it's very small - but non-zero to allow the gamma parameters
#to be fit).


#Convert abundance to biomass. First, clean benthic data to include only benthic organisms that emerged, and remove previous summaries. Then gather by taxon and add the dry mass data

ben_clean <- ben %>%
  select(-notes,-total,-daphnia,-isopodia,-hemiptera,-copepod,-snails,-oligochaetes) %>%
  gather(taxon,abundance,"ceratopogonidae":"chironomidae")

#plot ambient vs cage
plot_ben_cage_v_amb <- ben_clean %>% 
  mutate(date = mdy(date)) %>% 
  filter(trt2 != "Exclusion",
         date != "2017-06-02") %>% 
  ggplot(aes(x = reorder(taxon, -abundance), y = abundance, color = trt2, shape = trt2))+
  geom_point(size = 2,position = position_jitterdodge(dodge.width = 0.6,
                                                      jitter.width = 0),
             alpha = .8)+
  theme_bw()+
  facet_grid(date~., scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 10))+
  scale_y_log10() +
  scale_color_grey(start = 0.2, end = 0.7)+
  theme(legend.title = element_blank())+
  xlab("prey_taxon")+
  ylab("Number per benthic sample")+
  #coord_flip() +
  NULL


ggsave(plot_ben_cage_v_amb, file = "plot_ben_cage_v_amb.tiff", dpi = 600, width = 6.5, height =4.5)



#make data to analyze - total benthic mg
ben_dm_tot <- merge(ben_clean,ind_mg_mean) %>% 
  mutate(sample_m2 = 0.76*0.3,
         benthic_mg_dm_m2 = (abundance*mean_mgdm)/sample_m2,
         benthic_mg01 = 0.01 + benthic_mg_dm_m2)