library(brms)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(lubridate)
library(RCurl)
library(cowplot)
library(janitor)


# Load data and brms model (or re-run the model below) ---------------------------------------------------------------

diet_mgdm <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/diet_mgdm.csv"))
diet_brms <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/diet_brms.rds?raw=true"))

# Make response non-zero. Make dates look like dates. Fix prey names. Remove Quillback (has zero diet items) 
d <- diet_mgdm %>% 
  mutate(mgdm01 = 0.01 +mg_diet_dm,
         date2 = as.factor(date)) %>% 
  filter(species != "quillback") %>% 
  mutate(species = str_replace(species, "river_shiner","rivershiner")) %>% 
  mutate(species = fct_relevel(species,c("bluegill","spotfin"))) %>% 
  mutate(prey_taxon = str_replace_all(prey_taxon, c("unk_ins" = "unidentified",
                                                    "corixid" = "corixidae",
                                                    "terr_hemip"= "hemiptera",
                                                    "springtail"= "collembola",
                                                    "zoop" = "zooplankton",
                                                    "terr_hymen" = "hymenoptera",
                                                    "mayfly" = "ephemeroptera")))


d %>% 
  #filter(grepl("chiro",prey_taxon)) %>% 
  ggplot(aes(x = reorder(prey_taxon, -mg_diet_dm), y = number, color = method)) + 
  geom_point(size = 3, position = position_jitterdodge(dodge.width = .5,
                                                       jitter.width = 0),
             alpha = 0.5)+
  scale_y_log10()+
  facet_wrap(~species) +
  #coord_flip() +
  NULL


d %>% 
  #filter(number > 1) %>% 
  #filter(grepl("chiro",prey_taxon)) %>% 
  ggplot(aes(x = reorder(prey_taxon, -number), y = number + 0.1, color = method)) + 
  geom_point(size = 2,position = position_jitterdodge(dodge.width = 0.6,
                                                      jitter.width = 0),
             alpha = .9)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 10))+
  scale_y_log10()+
  scale_color_brewer(type = "qual")+
  facet_grid(species ~ .) +
  #coord_flip() +
  NULL



  
  