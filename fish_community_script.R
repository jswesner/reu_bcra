# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(brms)
library(ggridges)
library(RCurl)
library(ggrepel)
library(janitor)

# Load data ---------------------------------------------------------------

fish_community <-read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/fish_community.csv"))

#We estimated fish densities in the backwater on four dates, 
#once before the experiment began, and three times after. On 
#each date, we used a 4 ft by 20 ft seine, and collected fish 
#in two 20 ft hauls on each shore, both above and below the beaver 
#dam. A subset of fishes was kept for diet analysis via gastric 
#lavage. All other fishes were counted and released.


# Summarize fish collections ----------------------------------------------

fish_table <- fish_community %>% 
  group_by(species, date) %>%  
  summarize(total_abund = sum(abund),
            total_area = sum(area_sampled_m2)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  summarize(total_fish = sum(total_abund),
            total_area_m2 = sum(total_area)) %>% 
  mutate(prop_abund = round(total_fish/sum(total_fish),3)) %>% 
  arrange(-total_fish)

write.csv(fish_table, file = "fish_table.csv")


#plot fish community over time
fish_over_time <- fish_community %>% 
  group_by(species, date) %>%  
  summarize(total_abund = sum(abund),
            total_area = sum(area_sampled_m2)) %>% 
  filter(total_abund >0) %>% 
  ggplot(aes(x = reorder(species,total_abund), y = total_abund)) +
  geom_bar(stat = "identity")+
  facet_grid(.~date) +
  coord_flip() +
  xlab("Fish species")+
  ylab("Abundance")+
  theme_classic()

ggsave(fish_over_time, file = "fish_over_time.tiff", dpi = 600, width = 6.5, height = 2, units = "in")
