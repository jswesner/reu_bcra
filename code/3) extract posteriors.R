#load packages and data
source("src/packages.R")

# Extract conditional posteriors 
# diet model --------------------------------------------------------------

diet_brms <- readRDS("models/diet_brms.rds")
fish_totals <- read_csv(file = "data/raw_data/fish_totals.csv")

diet_brm_postpreds <- diet_brms$data %>%
  distinct(prey_stage, date2, species) %>% 
  add_epred_draws(diet_brms, dpar = T, re_formula = NULL) %>% 
  filter(.draw <= 1000) %>% 
  left_join(fish_totals %>% mutate(date2 = as.factor(date2))) %>% 
  mutate(population_epred = .epred*total_abund) %>% 
  arrange(-.epred) 

# taxon_stage <- diet_brm_postpreds %>% ungroup() %>% distinct(prey_taxon) %>% 
#   separate(prey_taxon, c("prey_taxon_only", "prey_stage"), remove = F)

diet_ind_postpreds_wide <- diet_brm_postpreds %>% 
  group_by(prey_stage, date2, species, .draw) %>%
  summarize(sum = sum(.epred)) %>% 
  pivot_wider(names_from = prey_stage, values_from = sum) %>%
  replace(is.na(.), 0) %>% 
  mutate(total = pa + not_pa,
         prop_pa = pa/total) %>% 
  mutate(data_level = "Per capita")

diet_pop_postpreds_wide <- diet_brm_postpreds %>% 
  group_by(prey_stage, date2, species, .draw) %>% 
  summarize(sum = sum(population_epred)) %>% 
  pivot_wider(names_from = prey_stage, values_from = sum) %>%
  replace(is.na(.), 0) %>% 
  mutate(total = pa + not_pa,
         prop_pa = pa/total) %>% 
  mutate(data_level = "Per population")

diet_community <- diet_pop_postpreds_wide %>% 
  group_by(date2, .draw) %>% 
  summarize(pa = sum(pa),
            total = sum(total)) %>% 
  mutate(prop_pa = pa/total) %>% 
  mutate(data_level = "Per community",
         species = "Community")

all_diet_posts <- bind_rows(diet_ind_postpreds_wide, diet_pop_postpreds_wide, diet_community) %>% 
  mutate(species = case_when(species == "spotfin" ~ "Spotfin Shiner",
                             species == "bluegill" ~ "Bluegill",
                             species == "largemouth" ~ "Largemouth Bass",
                             species == "rivershiner" ~ "River Shiner",
                             TRUE ~ species)) %>% 
  ungroup() %>% 
  mutate(data_level = fct_relevel(data_level, "Per capita", "Per population"))

saveRDS(all_diet_posts, file = "posteriors/all_diet_posts.rds")


# emergence model ---------------------------------------------------------
emerge_dm_model <- readRDS("models/emerge_dm_model.rds")

emerge_cond_posts <- emerge_reu_mg %>% 
  data_grid(date, trt2) %>% 
  add_epred_draws(emerge_dm_model, re_formula = NA) %>% 
  mutate(date = mdy(date),
         trt = case_when(trt2 == "ctrl" ~ "fish", TRUE ~ "no fish"))

saveRDS(emerge_cond_posts, file = "posteriors/emerge_cond_posts.rds")

# benthic model ---------------
brm_ben_m2 <- readRDS("models/brm_ben_m2.rds")

benthic_cond_posts <- ben_dm_tot %>% 
  data_grid(date, trt, taxon) %>% 
  add_epred_draws(brm_ben_m2, re_formula = NA) %>% 
  mutate(date = mdy(date),
         trt = case_when(trt == "ctrl" ~ "fish", TRUE ~ "no fish"))
  
saveRDS(benthic_cond_posts, file = "posteriors/benthic_cond_posts.rds")

# spider model ---------------
spiders_brm <- readRDS("models/spiders_brm.rds")

spiders_cond_posts <- spider_abund %>% 
  select(-trt) %>% 
  mutate(date = mdy(date)) %>%
  rename(trt = treatment) %>% 
  mutate(trt = case_when(trt == "fish" ~ "ctrl", TRUE ~ "exc")) %>% 
  data_grid(date, trt) %>% 
  add_epred_draws(spiders_brm, re_formula = NA) %>% 
  mutate(trt = case_when(trt == "ctrl" ~ "fish", TRUE ~ "no fish"))

saveRDS(spiders_cond_posts, file = "posteriors/spiders_cond_posts.rds")


