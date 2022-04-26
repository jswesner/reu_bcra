library(cheddar)
library(ggrepel)
library(tidyverse)


d_stages <- read_csv("data/derived_data/d.csv") 
d_no_stages <- read_csv("data/derived_data/d.csv") %>% 
  separate(prey_taxon, c("prey_taxon", "prey_stage")) 



# with stages -------------------------------------------------------------

diet_prop <- d_stages %>% 
  group_by(species) %>% 
  mutate(total_gdm = sum(mg_diet_dm)) %>% 
  group_by(species, prey_taxon, total_gdm) %>% 
  summarize(total_gdm_taxa = sum(mg_diet_dm)) %>% 
  mutate(proportion = total_gdm_taxa/total_gdm)


diet_prop_plot <- diet_prop %>% 
  ungroup() %>% 
  mutate(consumer = species,
         resource = prey_taxon,
         diet_post = proportion,
         resource = fct_relevel(as.factor(resource), c("chiro_l","coleo_a", "chiro_p"))) %>% 
  select(consumer, resource, diet_post) %>% 
  filter(resource != "zooplankton" & resource !="unidentified")

consumer_x <- diet_prop_plot %>% 
  select(consumer) %>% 
  distinct() %>% 
  drop_na() %>% 
  mutate(xend = rescale(as.numeric(as.factor(consumer)), to = c(0.25,0.75)))

resource_x <- diet_prop_plot %>% ungroup() %>% select(resource) %>% 
  distinct(resource) %>% 
  drop_na() %>% 
  mutate(trophic_prey = ifelse(resource == "fish_yoy", "vert_prey",
                               ifelse(resource == "seed", "plant_prey", "invertebrate_prey"))) %>%
  group_by(trophic_prey) %>% 
  drop_na() %>% 
  ungroup() %>% 
  distinct(resource) %>% 
  mutate(x = rescale(as.numeric(as.factor(resource)), to = c(0,1)),
         offset = rnorm(n(), 0,0.015),
         x1 = x + offset)

edges <- diet_prop_plot %>% 
  select(resource, consumer, diet_post) %>% 
  drop_na(resource) %>%
  left_join(consumer_x) %>% 
  left_join(resource_x) %>% 
  mutate(size = diet_post,
         xend1 = xend + offset,
         color = ifelse(grepl("chiro", resource),"chironomid","non-chironomid")) %>%
  mutate(diet_post_corrected = diet_post)   %>% 
  arrange(resource)  %>% 
  group_by(consumer) %>% 
  mutate(repeats = diet_post*100, 
         y = 1,
         yend = 2) %>% 
  ungroup() %>% 
  mutate(y = ifelse(resource=="fish_yoy", 1.5, 
                    ifelse(resource == "seed_unknown",0.5, y))) %>% 
  group_by(resource, consumer)  %>% 
  ungroup() 


nodes <- edges %>%
  ungroup() %>% 
  mutate(resource = str_replace(resource, "_","-")) %>% 
  unite(consumer_xy, c("consumer", "xend", "yend")) %>% 
  unite(resource_xy, c("resource", "x", "y")) %>% 
  select(resource_xy, consumer_xy) %>% 
  gather(key, value, "resource_xy":"consumer_xy") %>% 
  separate(value, c("taxon", "x","y"), sep = "_") %>% 
  separate(key, c("trophic","delete")) %>%
  select(-delete) %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         color = ifelse(grepl("chir", taxon),"chironomid","non-chironomid")) %>% 
  separate(taxon, c("taxon_ungrouped", "stage"), remove = F) %>% 
  replace_na(list(stage = "unknown"))

labels <- nodes %>% 
  distinct(taxon, .keep_all = T) %>%
  mutate(prey_label = ifelse(grepl("chiro",taxon), taxon,NA),
         label = ifelse(taxon == "chironomidae", "all chiro stages",
                        ifelse(taxon == "chiro-l", "larvae",
                               ifelse(taxon== "chiro-p", "pupae",
                                      ifelse(taxon == "chiro-a","adults",NA)))))




# make plot
library(ggthemes)

edges %>% 
  filter(diet_post_corrected != 0) %>% 
  arrange(consumer) %>% 
  ggplot() +
  geom_label_repel(data = labels, aes(x = x, y = y, label = label),
                   nudge_y = -.2, size = 3) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = diet_post_corrected,
                   alpha = diet_post_corrected)) +
  geom_point(data = nodes, aes(x = x, y = y, fill = stage),
             size=8, shape = 21) +
  #theme_classic() +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        text = element_text(size = 12)) +
  scale_fill_colorblind() +
  #scale_fill_viridis_d()+
  #scale_fill_brewer(type = "qual", palette = 4, direction = 1)+
  #scale_color_grey(start = 0.2, end = 0.8)+
  scale_y_continuous(breaks = c(0.5,1,1.5,2), labels = c("1","2","3","4"))+
  coord_cartesian(ylim = c(0.4,2.1)) +
  guides(size = guide_legend("Proportion of diet"),
         fill = guide_legend(""),
         color = F,
         alpha = F)+
  ylab("Trophic position")






# without stages ----------------------------------------------------------


diet_prop <- d_no_stages %>% 
  group_by(species) %>% 
  mutate(total_gdm = sum(mg_diet_dm)) %>% 
  group_by(species, prey_taxon, total_gdm) %>% 
  summarize(total_gdm_taxa = sum(mg_diet_dm)) %>% 
  mutate(proportion = total_gdm_taxa/total_gdm)


diet_prop_plot <- diet_prop %>% 
  ungroup() %>% 
  mutate(consumer = species,
         resource = prey_taxon,
         diet_post = proportion,
         resource = fct_relevel(as.factor(resource), c("chiro"))) %>% 
  select(consumer, resource, diet_post) %>% 
  filter(resource != "zooplankton" & resource !="unidentified")

consumer_x <- diet_prop_plot %>% 
  select(consumer) %>% 
  distinct() %>% 
  drop_na() %>% 
  mutate(xend = rescale(as.numeric(as.factor(consumer)), to = c(0.25,0.75)))

resource_x <- diet_prop_plot %>% ungroup() %>% select(resource) %>% 
  distinct(resource) %>% 
  drop_na() %>% 
  mutate(trophic_prey = ifelse(resource == "fish_yoy", "vert_prey",
                               ifelse(resource == "seed", "plant_prey", "invertebrate_prey"))) %>%
  group_by(trophic_prey) %>% 
  drop_na() %>% 
  ungroup() %>% 
  distinct(resource) %>% 
  mutate(x = rescale(as.numeric(as.factor(resource)), to = c(0,1)),
         offset = rnorm(n(), 0,0.015),
         x1 = x + offset)

edges <- diet_prop_plot %>% 
  select(resource, consumer, diet_post) %>% 
  drop_na(resource) %>%
  left_join(consumer_x) %>% 
  left_join(resource_x) %>% 
  mutate(size = diet_post,
         xend1 = xend + offset,
         color = ifelse(grepl("chiro", resource),"chironomid","non-chironomid")) %>%
  mutate(diet_post_corrected = diet_post)   %>% 
  arrange(resource)  %>% 
  group_by(consumer) %>% 
  mutate(repeats = diet_post*100, 
         y = 1,
         yend = 2) %>% 
  ungroup() %>% 
  mutate(y = ifelse(resource=="fish", 1.5, 
                    ifelse(resource == "seed",0.5, y))) %>% 
  group_by(resource, consumer)  %>% 
  ungroup() 


nodes <- edges %>%
  ungroup() %>% 
  mutate(resource = str_replace(resource, "_","-")) %>% 
  unite(consumer_xy, c("consumer", "xend", "yend")) %>% 
  unite(resource_xy, c("resource", "x", "y")) %>% 
  select(resource_xy, consumer_xy) %>% 
  gather(key, value, "resource_xy":"consumer_xy") %>% 
  separate(value, c("taxon", "x","y"), sep = "_") %>% 
  separate(key, c("trophic","delete")) %>%
  select(-delete) %>% 
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         color = ifelse(grepl("chir", taxon),"chironomid","non-chironomid"))

labels <- nodes %>% 
  distinct(taxon, .keep_all = T) %>%
  mutate(prey_label = ifelse(grepl("chiro",taxon), taxon,NA),
         label = ifelse(taxon == "chironomidae", "all chiro stages",
                        ifelse(taxon == "chiro-l", "larvae",
                               ifelse(taxon== "chiro-p", "pupae",
                                      ifelse(taxon == "chiro-a","adults",NA)))))


edges %>% 
  filter(diet_post_corrected != 0) %>% 
  arrange(consumer) %>% 
  ggplot() +
  geom_label_repel(data = labels, aes(x = x, y = y, label = label),
                   nudge_y = -.2, size = 3) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = diet_post_corrected,
                   alpha = diet_post_corrected)) +
  geom_point(data = nodes, aes(fill = color, x = x, y = y),
             size=8, shape = 21) +
  #theme_classic() +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        text = element_text(size = 12)) +
  #scale_fill_viridis_d()+
  #scale_fill_brewer(type = "qual", palette = 4, direction = 1)+
  #scale_color_grey(start = 0.2, end = 0.8)+
  scale_y_continuous(breaks = c(0.5,1,1.5,2), labels = c("1","2","3","4"))+
  coord_cartesian(ylim = c(0.4,2.1)) +
  guides(size = guide_legend("Proportion of diet"),
         fill = guide_legend(""),
         color = F,
         alpha = F)+
  ylab("Trophic position")

