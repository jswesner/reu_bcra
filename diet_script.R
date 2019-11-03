library(brms)
library(tidyverse)
library(ggridges)
library(lubridate)
library(RCurl)


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


# Plot raw data -----------------------------------------------------------

d %>% 
  ggplot(aes(x = date, y = mg_diet_dm + 0.01, color = prey_taxon)) +
  geom_point(position = position_jitter(), alpha =0.5) +
  facet_wrap(~species, scales = "free")+
  scale_y_log10()


# Bayesian model ----------------------------------------------------------

#diet_brms <- brm(mgdm01 ~ prey_taxon*date2*species, family = Gamma(link = "log"),
      #data = d,
      #prior = c(prior(normal(0,4), class = "Intercept"),
      #prior(normal(0,2), class = "b")),
      #cores = 4)
#saveRDS(diet_brms, file = "diet_brms.rds")

# Extract conditional posteriors ------------------------------------------

#make data frame to condition on
date2 <- unique(d$date2)
species <- unique(d$species)
prey_taxon <- unique(d$prey_taxon)

new_data <- expand_grid(species,date2, prey_taxon)
names <- new_data %>% unite(colnames,c(species, date2, prey_taxon)) #names for the columns of the fitted estimates below

#extract posteriors
fit_dir <- fitted(all_diet, newdata = new_data, re_formula = NA, summary = F)
colnames(fit_dir) <- names$colnames
as_tibble(fit_dir) 

# tidy posterior
post <- as_tibble(fit_dir) %>% 
  mutate(iter = 1:nrow(fit_dir)) %>% 
  gather(key, mg_dm_diet,-iter) %>% 
  separate(key, c("species","date2","prey_taxon"), sep = "_", extra = "merge")



# Plot posteriors ------------------------------------------

#collapse aggregate chiro stages to just one group
post_agg <- post %>% 
  pivot_wider(names_from = prey_taxon, 
              values_from = mg_dm_diet) %>% 
  mutate(chironomidae = chiro_l + chiro_a + chiro_p,
         coleoptera = coleo_ad + coleo_l) %>% 
  select(-chiro_a, -chiro_p, -chiro_l, -coleo_l, -coleo_ad) %>% 
  gather(prey_taxon, mg_dm_diet, c(-iter,-species,-date2)) %>% 
  group_by(species,date2,prey_taxon) %>% 
  ungroup()

#combine posts with correction
post1000_agg_0 <- post_agg %>% 
  ungroup() %>% 
  mutate(species = fct_relevel(species, c("bluegill","spotfin")))

#plot total mg prey
plot_mg_rawdata <- d %>% 
  select(id, species, date2, prey_taxon, mg_diet_dm) %>% 
  pivot_wider(names_from = prey_taxon,
              values_from = mg_diet_dm) %>% 
  mutate(chironomidae = chiro_l + chiro_p + chiro_a,
         coleoptera = coleo_ad + coleo_l) %>% 
  select(-chiro_l, -chiro_a, -chiro_p, -coleo_ad,-coleo_l) %>% 
  gather(prey_taxon, mg_dm_diet, c(-id,-species,-date2))
  
  
plot_mg_post <- post1000_agg_0 %>%
  filter(prey_taxon != "fish_yoy") %>% 
  mutate(mg_dm_diet_corrected = mg_dm_diet*correction)


plot_mg_diet <- 
  ggplot(data = plot_mg_post,aes(x = reorder(prey_taxon,mg_dm_diet), y = mg_dm_diet_corrected, 
             fill = date2,
             group = interaction(prey_taxon, date2))) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(-.8), size = 0.5,
               fatten = 1)+
  #scale_y_log10()+
  #ylim(c(0.1,30)) +
  facet_grid(.~species) +
  #geom_point(alpha = 0.5) +
  #geom_violin() + +
  #scale_y_log10()  +
  scale_fill_brewer(name = "Bayesian posterior") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top")+
  xlab("Prey taxon") +
  ylab("Fish diet (mg dry mass)")+
  geom_point(data = subset(plot_mg_rawdata,prey_taxon != "fish_yoy"),
             aes(x = reorder(prey_taxon,mg_dm_diet), y = mg_dm_diet,
                 group = interaction(prey_taxon, date2)),
             alpha = 0.4,
             size = 0.7,
             position = position_jitterdodge(dodge.width = -.8,
                                             jitter.width = 0,
                                             jitter.height = 0))

#ggsave(plot_mg_diet, file = "plot_mg_diet.tiff", dpi = 600, width = 7, height = 7, units = "in")



# Plot prop chiro and prop pupae 

#prop chiro prey
prop_chiro_with_fish <- as_tibble(post) %>% 
  pivot_wider(names_from = prey_taxon, 
              values_from = mg_dm_diet) %>% 
  mutate(chironomidae = chiro_l + chiro_a + chiro_p,
         coleoptera = coleo_ad + coleo_l) %>% 
  select(-chiro_a, -chiro_p, -chiro_l, -coleo_l, -coleo_ad) %>% 
  gather(prey_taxon, mg_dm_diet, c(-iter,-date2,-species)) %>% 
  mutate(chiro = ifelse(prey_taxon =="chironomidae", "chiro","non_chiro")) %>% 
  group_by(iter, species, date2, chiro) %>% 
  summarize(mg_dm_diet= sum(mg_dm_diet)) %>% 
  pivot_wider(names_from = chiro,
              values_from = mg_dm_diet) %>% 
  mutate(prop_chiro = chiro/(chiro + non_chiro),
         date_order = as.numeric(as.factor(date2)))


#prop chiro prey of non-fish prey
prop_chiro_without_fish <- as_tibble(post) %>% 
  filter(prey_taxon != "fish_yoy") %>% 
  pivot_wider(names_from = prey_taxon, 
              values_from = mg_dm_diet) %>% 
  mutate(chironomidae = chiro_l + chiro_a + chiro_p,
         coleoptera = coleo_ad + coleo_l) %>% 
  select(-chiro_a, -chiro_p, -chiro_l, -coleo_l, -coleo_ad) %>% 
  gather(prey_taxon, mg_dm_diet, c(-iter,-date2,-species)) %>% 
  mutate(chiro = ifelse(prey_taxon =="chironomidae", "chiro","non_chiro")) %>% 
  group_by(iter, species, date2, chiro) %>% 
  summarize(mg_dm_diet= sum(mg_dm_diet)) %>% 
  pivot_wider(names_from = chiro,
              values_from = mg_dm_diet) %>% 
  mutate(prop_chiro = chiro/(chiro + non_chiro),
         date_order = as.numeric(as.factor(date2))) %>% 
  ungroup() %>% 
  mutate(species = str_replace_all(species, c("bluegill" = "Bluegill",
                                              "largemouth" = "Largemouth Bass",
                                              "rivershiner" = "River Shiner",
                                              "spotfin" = "Spotfin Shiner")),
         species = fct_relevel(species, c("Bluegill", "Spotfin Shiner"))) 

#prop chiro that are pupae
prop_chiro_pupae <- as_tibble(post) %>% 
  filter(prey_taxon == "chiro_l" | prey_taxon == "chiro_a" | prey_taxon == "chiro_p") %>% 
  pivot_wider(names_from = prey_taxon, 
              values_from = mg_dm_diet) %>% 
  mutate(chironomidae = chiro_l + chiro_a + chiro_p) %>% 
  mutate(prop_pa = 1-(chiro_l/chironomidae)) %>% 
  mutate(species = str_replace_all(species, c("bluegill" = "Bluegill",
                                              "largemouth" = "Largemouth Bass",
                                              "rivershiner" = "River Shiner",
                                              "spotfin" = "Spotfin Shiner")),
         species = fct_relevel(species, c("Bluegill", "Spotfin Shiner")))


#combine and plot in two panels
both <- prop_chiro_pupae %>% left_join(prop_chiro_without_fish) %>% 
  select(iter, species, date2, prop_pa, prop_chiro) %>% 
  gather(facet, proportion, c("prop_pa", "prop_chiro")) %>% 
  mutate(facet = str_replace_all(facet, c("prop_pa" = "Proportion pupa + adult",
                                          "prop_chiro" = "Proportion chironomid")))


plot_prop_chiro_pup <- both %>%  
  filter(species == "Bluegill" | species == "Spotfin Shiner") %>% 
  ggplot(aes(x = proportion, y = fct_rev(date2), fill = species))+
  geom_density_ridges() +
  scale_fill_grey()+
  theme_ridges() +
  theme(legend.title = element_blank())+
  ylab("Sample Date")+
  xlab("Proportion in diets")+
  facet_wrap(~facet)


#ggsave(plot_prop_chiro_pup, file = "plot_prop_chiro_pup.tiff", dpi = 600, width = 8, height = 5, units = "in")



# Summarize posteriors ----------------------------------------------------
#summary stats of dry mass in fish diets
diet_post <- post_agg %>% 
  group_by(species, date2, prey_taxon) %>% 
  summarize(mean = mean(mg_dm_diet),
            median = median(mg_dm_diet),
            sd = sd(mg_dm_diet),
            low95 = quantile(mg_dm_diet, probs = 0.025),
            high95 = quantile(mg_dm_diet, probs = 0.975)) %>% 
  arrange(-mean) 

#summary stats by proportion of diet items in fish diets

#no stages
sum_diet <- post_agg %>% 
  group_by(species, date2, iter) %>%
  summarize(sum = sum(mg_dm_diet))

diet_post_proportion_nostage <- post_agg %>% 
  left_join(sum_diet) %>% 
  group_by(species, date2, prey_taxon) %>% 
  mutate(prop = mg_dm_diet/sum) %>% 
  summarize(mean = mean(prop),
            median = median(prop),
            var = var(prop), 
            precision = 1/var,
            sd = sd(prop),
            low95 = quantile(prop, probs = 0.025),
            high95 = quantile(prop, probs = 0.975)) %>% 
  arrange(-mean) %>% 
  mutate(aggregation = "no_stage_structure")

#with stages
sum_diet_stage <- post %>% 
  group_by(species, date2, iter) %>%
  summarize(sum = sum(mg_dm_diet))

diet_post_proportion_stage <- post %>% 
  left_join(sum_diet_stage) %>% 
  group_by(species, date2, prey_taxon) %>% 
  mutate(prop = mg_dm_diet/sum) %>% 
  summarize(mean = mean(prop),
            median = median(prop),
            var = var(prop), 
            precision = 1/var,
            sd = sd(prop),
            low95 = quantile(prop, probs = 0.025),
            high95 = quantile(prop, probs = 0.975)) %>% 
  arrange(-mean) %>%
  mutate(aggregation = "stage_structure")

#combine
diet_post_proportion <- rbind(diet_post_proportion_nostage, diet_post_proportion_stage)
saveRDS(diet_post_proportion, file = "diet_post_proportion.rds")

#summary stats of prop chiro
#WITH FISH_yoy
prop_chiro_with_fish %>% 
  group_by(species, date2) %>% 
  summarize(mean = mean(prop_chiro),
            median = median(prop_chiro),
            sd = sd(prop_chiro),
            low95 = quantile(prop_chiro, probs = 0.025),
            high95 = quantile(prop_chiro, probs = 0.975)) %>% 
  arrange(-mean) 

#WITHOUT FISH_yoy
prop_chiro_without_fish %>% 
  group_by(species, date2) %>% 
  summarize(mean = mean(prop_chiro),
            median = median(prop_chiro),
            sd = sd(prop_chiro),
            low95 = quantile(prop_chiro, probs = 0.025),
            high95 = quantile(prop_chiro, probs = 0.975)) %>% 
  arrange(-mean) 





#summary stats of prop chiro
#WITH FISH_yoy
prop_chiro_pupae %>% 
  group_by(species, date2) %>% 
  summarize(mean = mean(prop_pa),
            median = median(prop_pa),
            sd = sd(prop_pa),
            low95 = quantile(prop_pa, probs = 0.025),
            high95 = quantile(prop_pa, probs = 0.975)) %>% 
  arrange(-mean)

