library(brms)
library(tidyverse)
library(ggridges)
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


# Plot raw data -----------------------------------------------------------

d %>% 
  ggplot(aes(x = date, y = mg_diet_dm + 0.01, color = prey_taxon)) +
  geom_point(position = position_jitter(), alpha =0.5) +
  facet_wrap(~species, scales = "free")+
  scale_y_log10()

#compare methods
plot_diet_method <- d %>% 
  #filter(number > 1) %>% 
  #filter(grepl("chiro",prey_taxon)) %>% 
  ggplot(aes(x = reorder(prey_taxon, -mg_diet_dm), y = mg_diet_dm + 0.1, color = method,
             shape = method)) + 
  geom_point(size = 2,position = position_jitterdodge(dodge.width = 0.6,
                                                      jitter.width = 0),
             alpha = .9)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 10))+
  scale_y_log10()+
  scale_color_grey(start = 0.2, end = 0.7)+
  facet_grid(species ~ .) +
  xlab("prey_taxon")+
  ylab("Dry mass per stomach (mg)")+
  #coord_flip() +
  NULL

ggsave(plot_diet_method, file = "plot_diet_method.tiff", dpi = 600, width = 6.5, height = 7)

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

#make a correction column that contains zeros and ones to indicate 
#whether a prey item was present

correction <- d %>% 
  select(species, date, prey_taxon, mg_diet_dm) %>%
  group_by(species, date, prey_taxon) %>% 
  summarize(sum_dm = sum(mg_diet_dm)) %>% 
  pivot_wider(names_from = prey_taxon, 
              values_from = sum_dm) %>% 
  mutate(chironomidae = chiro_l + chiro_a + chiro_p,
         coleoptera = coleo_ad + coleo_l) %>% 
  select(-chiro_a, -chiro_p, -chiro_l, -coleo_l, -coleo_ad) %>% 
  gather(prey_taxon, sum_dm, c(-species,-date)) %>% 
  mutate(correction = ifelse(sum_dm >0, 1, 0)) %>% 
  select(-sum_dm)
  

#combine posts with correction
post1000_agg_0 <- post_agg %>% 
  left_join(correction) %>%
  ungroup() %>% 
  mutate(species = fct_relevel(species, c("bluegill","spotfin")))

#total mg prey raw data
plot_mg_rawdata <- d %>% 
  select(id, species, date2, prey_taxon, mg_diet_dm) %>% 
  pivot_wider(names_from = prey_taxon,
              values_from = mg_diet_dm) %>% 
  mutate(chironomidae = chiro_l + chiro_p + chiro_a,
         coleoptera = coleo_ad + coleo_l) %>% 
  select(-chiro_l, -chiro_a, -chiro_p, -coleo_ad,-coleo_l) %>% 
  gather(prey_taxon, mg_dm_diet, c(-id,-species,-date2))
  
#posterior plotting data  excluding fish_yoy
plot_mg_post <- post1000_agg_0 %>%
  filter(prey_taxon != "fish_yoy") %>% 
  mutate(mg_dm_diet_corrected = mg_dm_diet*correction)

#plot_without fish_yoy
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
                                             jitter.height = 0))+
  guides(fill = F)

#posterior_with_fishyoy only
plot_mg_post_withfish <- post1000_agg_0 %>%
  filter(prey_taxon == "fish_yoy") %>% 
  mutate(mg_dm_diet_corrected = mg_dm_diet*correction)

#plot_with_fish_yoy only
plot_mg_diet_fishyoy <- 
  ggplot(data = plot_mg_post_withfish,aes(x = prey_taxon, y = mg_dm_diet, 
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
  coord_flip(ylim = c(0,1100)) +
  theme_classic() +
  theme(legend.position = "top",
        axis.title.x = element_blank())+
  xlab("") +
  geom_point(data = subset(plot_mg_rawdata,prey_taxon == "fish_yoy"),
             aes(x = reorder(prey_taxon,mg_dm_diet), y = mg_dm_diet,
                 group = interaction(prey_taxon, date2)),
             alpha = 0.4,
             size = 0.7,
             position = position_jitterdodge(dodge.width = -.8,
                                             jitter.width = .3,
                                             jitter.height = 0))

#combine plots with and without fish_yoy so that scales for insects are not overwhelmed by
#fish yoy

plot_all_diet <- plot_grid(plot_mg_diet_fishyoy, plot_mg_diet,
                           rel_heights = c(0.3,1),
                           align = "v",
                           ncol = 1)



ggsave(plot_all_diet , file = "plot_mg_diet.tiff", dpi = 600, width = 7, height = 10, units = "in")



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
  mutate(facet = str_replace_all(facet, c("prop_pa" = "b) Proportion pupa + adult",
                                          "prop_chiro" = "a) Proportion chironomid")))


plot_prop_chiro_pup <- both %>%  
  #filter(species == "Bluegill" | species == "Spotfin Shiner") %>% 
  ggplot(aes(x = proportion, y = fct_rev(date2), fill = fct_relevel(species,"River Shiner", "Largemouth Bass")))+
  geom_density_ridges(alpha = 1) +
  #scale_fill_grey(start = 1, end = 0)+
  scale_fill_manual(values = c("grey100","grey92","grey50","grey0"))+
  theme_classic() +
  theme(legend.title = element_blank(),
        strip.text = element_text(angle = 0, hjust = 0))+
  ylab("Sample Date")+
  xlab("Proportion in diets")+
  facet_wrap(~facet)


ggsave(plot_prop_chiro_pup, file = "plot_prop_chiro_pup.tiff", dpi = 600, width = 8, height = 5, units = "in")



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
  group_by(species,date2, prey_taxon) %>% 
  mutate(prop = mg_dm_diet/sum) %>% 
  summarize(mean = mean(prop),
            median = median(prop),
            var = var(prop), 
            precision = 1/var,
            sd = sd(prop),
            low95 = quantile(prop, probs = 0.025),
            high95 = quantile(prop, probs = 0.975)) %>% 
  arrange(species,-mean) %>% 
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
  arrange(species,date2,-mean) 


#summary stats of prop chiro
#WITH FISH_yoy
prop_chiro_pupae_tbl <- prop_chiro_pupae %>% 
  group_by(species, date2) %>% 
  summarize(mean = mean(prop_pa),
            median = median(prop_pa),
            sd = sd(prop_pa),
            low95 = quantile(prop_pa, probs = 0.025),
            high95 = quantile(prop_pa, probs = 0.975)) %>% 
  mutate_if(is.numeric, round, digits = 2) %>% 
  arrange(species,date2,-mean) 

write.csv(prop_chiro_pupae_tbl, file = "prop_chiro_pupae.csv")

#probability of a difference between species
prop_chiro_pupae %>% 
  select(iter,species, date2, prop_pa) %>% 
  pivot_wider(names_from = species,
              values_from = prop_pa) %>% 
  clean_names() %>% 
  group_by(date2) %>% 
  mutate(diff_spot_blue = spotfin_shiner - bluegill) %>% 
  summarize(prob_spot_greaterthan_blue = sum(diff_spot_blue>0)/length(diff_spot_blue))

#magnitude of difference
prop_chiro_pupae %>% 
  select(iter,species, date2, prop_pa) %>% 
  pivot_wider(names_from = species,
              values_from = prop_pa) %>% 
  clean_names() %>% 
  group_by(date2) %>% 
  mutate(diff_spot_blue = spotfin_shiner - bluegill) %>% 
  summarize(mean = mean(diff_spot_blue),
            low95 = quantile(diff_spot_blue, probs = 0.025),
            upper95 = quantile(diff_spot_blue, probs = 0.975))

#overall mean prop pupae
diet_mgdm %>% 
  select(species, date, prey_taxon, mg_diet_dm) %>% 
  filter(grepl("chiro",prey_taxon)) %>% 
  summarize(tot_chiro_raw = sum(mg_diet_dm))

prop_chiro_pupae %>% 
  summarize(mean = mean(prop_pa),
            low95 = quantile(prop_pa, probs = 0.025),
            upper95 = quantile(prop_pa, probs = 0.975))

# proportion of pupae + adults
both %>% 
  group_by(species, date2, facet) %>% 
  summarize(mean = mean(proportion),
            median = median(proportion),
            sd = sd(proportion),
            low95 = quantile(proportion, probs = 0.025),
            hihg95 = quantile(proportion, probs = 0.975)) %>% 
  arrange(desc(facet))



# Summarize raw data ------------------------------------------------------

#freq_occurence
diet_mgdm %>% 
  mutate(chiro_ind = ifelse(grepl("chiro", prey_taxon) & mg_diet_dm > 0, 1, 0)) %>% 
  group_by(id) %>% 
  summarize(chiro_present = sum(chiro_ind)) %>% 
  mutate(chiro_pres = ifelse(chiro_present > 0, 1,0)) 
  summarize(prop_pres = sum(chiro_pres)/length(id))


sp_id <- diet_mgdm %>% select(id, species) %>% distinct(id, .keep_all = T)

freq_occur <- diet_mgdm %>% 
  filter(species != "quillback") %>% 
  select(mg_diet_dm, id, prey_taxon) %>% 
  pivot_wider(names_from = prey_taxon,
              values_from = mg_diet_dm) %>% 
  mutate(chiro = chiro_l + chiro_a + chiro_p) %>% 
  select(-chiro_l, -chiro_a, -chiro_p) %>% 
  gather(prey_taxon, mg_diet_dm, c(-id)) %>% 
   # mutate(chiro_ind = ifelse(grepl("chiro", prey_taxon) & mg_diet_dm > 0, 1, 0)) %>% 
  mutate(taxon_pres = ifelse(mg_diet_dm > 0, 1,0))%>% 
  left_join(sp_id) %>% 
  group_by(prey_taxon, species) %>% 
  summarize(prop_pres = sum(taxon_pres)/length(id)) %>% 
  arrange(-prop_pres, species) %>% 
  pivot_wider(names_from = species, 
              values_from = prop_pres) %>% 
  select(prey_taxon, bluegill, spotfin, largemouth, river_shiner)

write.csv(freq_occur, file = "freq_occur.csv")  

diet_mgdm %>% 
  distinct(id, .keep_all = T) %>% 
  select(species, date) %>% 
  group_by(species, date) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = date,
              values_from = n)
