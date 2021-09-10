#load packages and data
source("code/1) load data.R")

all_diet_posts <- readRDS(file = "posteriors/all_diet_posts.rds")


# Summarize diets ---------------------------------------------------

# species averages
species_global <- all_diet_posts %>% 
  pivot_longer(cols = c(-date2, -species, -.draw, -data_level)) %>%  
  group_by(data_level, .draw, name) %>% 
  summarize(value = mean(value, na.rm = T)) %>% 
  mutate(species = "All Species") %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  bind_rows(all_diet_posts %>% 
              pivot_longer(cols = c(-date2, -species, -.draw, -data_level)) %>%  
              group_by(data_level, species, .draw, name) %>% 
              summarize(value = mean(value, na.rm = T)) %>% 
              pivot_wider(names_from = name, values_from = value)) %>% 
  pivot_longer(cols = c(-species, -.draw, -data_level)) 

species_means_global <- species_global %>% 
  group_by(data_level, species, name)  %>% 
  mutate(value = case_when(name == "total" ~ value*1000, TRUE ~ value)) %>% 
  summarize(median = median(value),
            mean = mean(value),
            sd = sd(value),
            lower = quantile(value, probs = 0.025, na.rm = T),
            upper = quantile(value, probs = 0.975, na.rm = T)) %>% 
  mutate_if(is.numeric, ~round(.,2)) %>% 
  filter(name %in% c("prop_pa", "total"))

write.csv(species_means_global, file = "tables/species_means_global.csv", row.names = F)

# species by date averages
species_by_date <- all_diet_posts %>% 
  pivot_longer(cols = c(-date2, -species, -.draw, -data_level)) %>%  
  group_by(data_level, .draw, date2, name) %>% 
  summarize(value = mean(value, na.rm = T)) %>% 
  mutate(species = "All Species") %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  bind_rows(all_diet_posts) %>% 
  pivot_longer(cols = c(-date2, -species, -.draw, -data_level)) 

species_means_bydate <- species_by_date %>% 
  group_by(data_level, date2, species, name)  %>% 
  mutate(value = case_when(name == "total" ~ value*1000, TRUE ~ value)) %>% 
  summarize(median = median(value),
            mean = mean(value),
            sd = sd(value),
            lower = quantile(value, probs = 0.025, na.rm = T),
            upper = quantile(value, probs = 0.975, na.rm = T)) %>% 
  mutate_if(is.numeric, ~round(.,2))  %>% 
  filter(name %in% c("prop_pa", "total"))

write.csv(species_means_bydate, file = "tables/species_means_bydate.csv", row.names = F)
