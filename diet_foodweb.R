library(tidyverse)
library(ggrepel)
library(cheddar)
library(scales)
library(RCurl)
library(viridis)

# Load data ---------------------------------------------------------------
#raw diet data
diet_mgdm <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/diet_mgdm.csv"))
#posterior of diet proportions
diet_brms <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/diet_post_proportion.rds?raw=true")) %>% 
  mutate(aggregation = ifelse(grepl("no_stage", aggregation), "a) no_stage_structure", "b) stage_structure"))


# Load functions ----------------------------------------------------------


# Create data to plot -----------------------------------------------------

# grid that assigns 0 or 1 to a predator-prey interaction based on raw data (use to remove links that didn't occur in samples)
# This is necessary because the bayesian model can only have values >0 and <1, but not true zeros.
chiro_coleo_agg <- diet_mgdm %>% 
  ungroup() %>% 
  mutate(consumer = species,
         resource = prey_taxon) %>% 
  select(consumer, date, resource, number, id) %>% 
  group_by(consumer, date, resource) %>% 
  mutate(sum = sum(number)) %>%
  ungroup() %>% 
  select(-number) %>% 
  group_by(consumer, date, resource) %>% 
  summarize(sum = sum(sum)) %>% 
  pivot_wider(names_from = resource,
              values_from = sum)  %>% 
  mutate(coleoptera = coleo_ad + coleo_l,
         chironomidae = chiro_a + chiro_l + chiro_p) %>% 
  select(consumer, date, coleoptera, chironomidae) %>% 
  gather(resource, number, "coleoptera":"chironomidae")



d_0 <- diet_mgdm %>% 
  ungroup() %>% 
  mutate(consumer = as.character(species),
         resource = as.character(prey_taxon), 
         date = as.character(date)) %>% 
  select(consumer, date, resource, number, id) %>% 
  group_by(consumer, date, resource) %>% 
  summarize(number = sum(number)) %>%
  ungroup() %>% 
  bind_rows(chiro_coleo_agg) %>% 
  mutate(correction = ifelse(number == 0, 0, 1),
         consumer = str_replace(consumer, "river_shiner","rivershiner")) %>% 
  select(-number)
  

#combine all food webs
all <- diet_brms %>% 
  ungroup() %>% 
  mutate(consumer = species,
         resource = prey_taxon,
         date = date2,
         diet_post = mean,
         resource = fct_relevel(as.factor(resource), c("chiro_l","chiro_p"))) %>% 
  select(consumer, resource, date, diet_post, aggregation, sd) %>% 
  filter(resource != "zooplankton" &resource !="unidentified")


#assign x positions to consumers. Make them evenly spaced between 0.25 and 0.75 (so the top is narrower)
consumer_x <- all %>% select(consumer) %>% 
  distinct() %>% 
  drop_na() %>% 
  mutate(xend = rescale(as.numeric(as.factor(consumer)), to = c(0.25,0.75))) 

#assign x positions to consumers. Make them evenly spaced between 0 and 1 
resource_x <- all %>% ungroup() %>% select(resource, date, aggregation) %>% 
  group_by(date, aggregation) %>% 
  distinct(resource) %>% 
  drop_na() %>% 
  mutate(trophic_prey = ifelse(resource == "fish_yoy", "vert_prey",
                               ifelse(resource == "seed", "plant_prey", "invertebrate_prey"))) %>%
  group_by(trophic_prey, aggregation) %>% 
  drop_na() %>% 
  ungroup() %>% 
  distinct(resource) %>% 
  mutate(resource = fct_relevel(as.factor(resource),"chironomidae","chiro_l","collembola","copepod","corixidae",
                                "ephemeroptera", "ostracod","fish_yoy","chiro_p","hymenoptera","seed",
                                "sciaridae","hemiptera","coleo_ad","coleo_l","spider", "chiro_a"))  %>%
  mutate(x = rescale(as.numeric(as.factor(resource)), to = c(0,1)),
         offset = rnorm(n(), 0,0.015),
         x1 = x + offset)
  


#no uncertainty
edges <- all %>% 
  select(resource, consumer, date, diet_post, aggregation, sd) %>% 
  drop_na(resource) %>%
  left_join(consumer_x) %>% 
  left_join(resource_x) %>% 
  mutate(size = diet_post,
         xend1 = xend + offset,
         color = ifelse(grepl("chiro", resource),"chironomid","non-chironomid")) %>%  
  group_by(aggregation, date) %>%
  left_join(d_0) %>% 
  mutate(diet_post_corrected = diet_post*correction)   %>% 
  arrange(resource)  %>% 
  group_by(consumer) %>% 
  mutate(repeats = diet_post*100, 
         y = 1,
         yend = 2) %>% 
  ungroup() %>% 
  mutate(y = ifelse(resource=="fish_yoy", 1.5, 
                    ifelse(resource == "seed",0.5, y))) %>% 
  group_by(resource, consumer, date, aggregation)  %>% 
  ungroup() 






nodes <- edges %>%
  ungroup() %>% 
  mutate(resource = str_replace(resource, "_","-")) %>% 
  unite(consumer_xy, c("consumer", "xend", "yend")) %>% 
  unite(resource_xy, c("resource", "x", "y")) %>% 
  select(resource_xy, consumer_xy, date, aggregation) %>% 
  gather(key, value, "resource_xy":"consumer_xy") %>% 
  separate(value, c("taxon", "x","y"), sep = "_") %>% 
  separate(key, c("trophic","delete")) %>%
  select(-delete) %>% 
  drop_na() %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         color = ifelse(grepl("chir", taxon),"chironomid","non-chironomid"))

labels <- nodes %>% 
  group_by(date, aggregation) %>% 
  distinct(taxon, .keep_all = T) %>%
  mutate(prey_label = ifelse(grepl("chiro",taxon), taxon,NA),
         label = ifelse(taxon == "chironomidae", "all chiro stages",
               ifelse(taxon == "chiro-l", "larvae",
                      ifelse(taxon== "chiro-p", "pupae",
                             ifelse(taxon == "chiro-a","adults", NA))))) %>%
  mutate(taxon_short = str_sub(taxon, 1,3),
         taxon_short = case_when(taxon == "chiro-p" ~ "chi-p",
                                 taxon == "chiro-a" ~ "chi-a",
                                 taxon == "chiro-l" ~ "chi-l",
                                 taxon == "coleo-l" ~ "col-l",
                                 taxon == "coleo-ad" ~ "col-a",
                                 taxon == "collembola" ~ "coll", 
                                 TRUE ~ taxon_short),
         labels_fish = ifelse(taxon == "bluegill", "BG",
                ifelse(taxon == "largemouth","LM",
                       ifelse(taxon == "rivershiner","RS",
                              ifelse(taxon == "spotfin", "SP",NA)))))


plot_foodweb <- edges %>% 
  arrange(consumer) %>% 
ggplot() +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = diet_post_corrected, 
                   alpha = correction))  +
  facet_grid(date ~ aggregation, scales = "free")+
  geom_point(data = nodes, aes(x = x, y = y, fill = color, shape = color),
             size=8)+
  geom_label_repel(data = labels, aes(x = x, y = y, label = label),
                   nudge_y = -.3, size = 6) +
  geom_label(data = labels, aes(x = x, y = y, label = labels_fish),
                   nudge_y = .2,
                   nudge_x = 0, size = 6) +
  guides(size = guide_legend("Proportion of diet"),
         color = guide_legend("") ,
         alpha = F) +
  #theme_classic() +
  theme(panel.background = element_rect(color = "black", fill = "white"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), 
    text = element_text(size = 25)) +
  scale_shape_manual(values = c(22,21))+
  scale_fill_grey() +
  #scale_color_brewer(type = "qual", palette = 4, direction = 1)+
  #scale_color_grey(start = 0.2, end = 0.8)+
  scale_y_continuous(breaks = c(0.5,1,1.5,2), labels = c("1","2","3","4"))+
  coord_cartesian(ylim = c(0.4,2.5)) +
  ylab("Trophic position")

plot_foodweb

ggsave(plot_foodweb, file = "plot_foodweb.tiff", dpi = 600, width = 14, height = 16, units = "in")
#



plot_foodweb_key <- edges %>% 
  arrange(consumer) %>% 
  filter(date == "2017-06-13") %>% 
  ggplot()  +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
               color = "grey80", size = 0.3) +
  facet_grid(.~aggregation, scales = "free")+
  geom_point(data = nodes %>% filter(date == "2017-06-13"), aes(x = x, y = y, fill = color, shape = color),
             size=4)+
  geom_label_repel(data = labels %>% filter(date == "2017-06-13"), aes(x = x, y = y, label = taxon_short),
                   size = 3) +
  guides(size = guide_legend("Proportion of diet"),
         color = guide_legend("") ,
         alpha = F) +
  #theme_classic() +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_shape_manual(values = c(22,21))+
  scale_fill_grey() +
  #scale_color_brewer(type = "qual", palette = 4, direction = 1)+
  #scale_color_grey(start = 0.2, end = 0.8)+
  scale_y_continuous(breaks = c(0.5,1,1.5,2), labels = c("1","2","3","4"))+
  coord_cartesian(ylim = c(0.4,2.1)) +
  ylab("Trophic position") +
  guides(color = F,
         shape = F,
         fill = F)

plot_foodweb_key
ggsave(plot_foodweb_key, file = "plot_foodweb_key.jpg", dpi = 500, width = 6, height = 3)

table_taxon_key <- labels %>% 
  ungroup() %>% 
  select(taxon, taxon_short) %>% 
  distinct(taxon, .keep_all = T) %>% 
  print(n = Inf)

write.csv(table_taxon_key, file = "table_taxon_key.csv",
          row.names = F)


# 
# 
# #Simulate uncertainty - needs work
# edges_sim <- edges %>%  
#   mutate(repeats = diet_post*300) %>% 
#   ungroup() %>% 
#   group_by(resource, consumer, date, aggregation) %>% 
#   slice(rep(1:n(), each = repeats)) %>% 
#   mutate(offset = rnorm(n(), 0,0.0005),
#          x1 = x+offset,
#          xend1 = xend + offset) %>% 
#   mutate(offset = rnorm(n(), 0,0.009),
#          x1 = x+offset,
#          xend1 = xend + offset,
#          size = diet_post)
# 
# 
# 
# 






