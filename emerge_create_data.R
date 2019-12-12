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

#final dataset to use with Bayesian emergence model
emerge_reu_mg_long <- emerge_reu_mg %>% 
  select(date,id, trt2, days, area, 
         chiro_mg_dm, doli_mg_dm, 
         tric_mg_dm, odo_mg_dm, ephem_mg_dm) %>% 
  mutate(trt3 = ifelse(date == "05/28/17", "ctrl", trt2)) %>% 
  gather(taxon, mg_dm, c(-date, -id, -trt2, -days, -area))

