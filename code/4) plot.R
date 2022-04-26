#load packages
source("src/packages.R")

#load posteriors
all_diet_posts <- readRDS(file = "posteriors/all_diet_posts.rds") %>% 
  mutate(data_level_ordered = case_when(data_level == "Per community" ~ "a) Per Community",
                                        data_level == "Per population" ~ "b) Per Population",
                                        TRUE ~ "c) Per capita"))
emerge_cond_posts <- readRDS(file = "posteriors/emerge_cond_posts.rds")
benthic_cond_posts <- readRDS("posteriors/benthic_cond_posts.rds")
spiders_cond_posts <- readRDS("posteriors/spiders_cond_posts.rds")

# load data
d <- read_csv("data/derived_data/d.csv")                               # diet data by taxa
d_tot <- read_csv("data/derived_data/d_tot.csv")                       # diet data analyzed
emerge_reu_mg <- read_csv("data/derived_data/emerge_reu_mg.csv") %>%   # emergence data
  mutate(date = mdy(date),
         trt2 = case_when(trt2 == "ctrl" ~ "fish", TRUE ~ "no fish"))
ben_dm_tot <- read_csv("data/derived_data/ben_dm_tot.csv")             # emergence data
spider_abund <- read_csv("data/derived_data/spider_abund.csv")  %>%    # spider data
  mutate(trt2 = case_when(trt == "ctrl" ~ "fish", TRUE ~ "no fish"))

# Plot raw data -----------------------------------------------------------

d %>% 
  ggplot(aes(x = date, y = mg_diet_dm, color = method)) +
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

ggsave(plot_diet_method, file = "plots/plot_diet_method.tiff", dpi = 600, width = 6.5, height = 7)


# compare ambient versus emergence cages
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

ggsave(plot_emerge_cage_v_amb, file = "plots/plot_emerge_cage_v_amb.tiff", dpi = 600, width = 6.5, height = 7)




# Plot diet --------------------------
prop_plot <- all_diet_posts %>% 
  filter(.draw <= 300) %>% 
  ungroup() %>% 
  mutate(species = fct_relevel(species, "Community", "Spotfin Shiner", "Bluegill"),
         date = ymd(date2),
         data_level_ordered = case_when(grepl("a)", data_level_ordered) ~ "d) Per Community",
                                        grepl("b)", data_level_ordered) ~ "e) Per Population",
                                        TRUE ~ "f) Per capita"),
         dodge_manual = case_when(species == "Spotfin Shiner" ~ 0.5,
                                  species == "Bluegill" ~ 0.25,
                                  species == "Community" ~ 0,
                                  species == "Largemouth Bass" ~ -0.25,
                                  TRUE ~ 0)) %>% 
  ggplot(aes(x = date + dodge_manual, y = prop_pa, fill = species)) +
  theme_tidybayes() +
  geom_line(aes(group = interaction(species, .draw), color = species), alpha = 0.01) +
  stat_pointinterval(aes(group = interaction(species, date), color = species), 
               alpha = 0.5) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  scale_y_log10(labels = comma, limits = c(0.001, 1)) +
  facet_wrap(~data_level_ordered) +
  coord_cartesian(xlim = c(as.Date("2017-06-12"), as.Date("2017-06-30"))) +
  labs(y = "Proportion pupae + adults",
       x = "", 
       color = "Species",
       fill = "Species") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  theme(strip.text.x = element_text(hjust = -0.01))

total_plot <- all_diet_posts %>% 
  filter(.draw <= 300) %>% 
  ungroup() %>% 
  mutate(species = fct_relevel(species, "Community", "Spotfin Shiner", "Bluegill"),
         data_level = fct_relevel(data_level, "Per community", "Per population"),
         date = ymd(date2),
         dodge_manual = case_when(species == "Spotfin Shiner" ~ 0.5,
                                  species == "Bluegill" ~ 0.25,
                                  species == "Community" ~ 0,
                                  species == "Largemouth Bass" ~ -0.25,
                                  TRUE ~ 0)) %>% 
  ggplot(aes(x = date + dodge_manual, y = total, fill = species)) +
  theme_tidybayes() +
  geom_line(aes(group = interaction(species, .draw), color = species), alpha = 0.01) +
  stat_pointinterval(aes(group = interaction(species, date), color = species), 
                     alpha = 0.5) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  facet_wrap(~data_level_ordered) +
  scale_y_log10(labels = comma) +
  labs(y = "Total prey eaten (mgDM)",
       x = "",
       color = "Species",
       fill = "Species") +
  coord_cartesian(xlim = c(as.Date("2017-06-12"), as.Date("2017-06-30"))) + 
  theme(strip.text.x = element_text(hjust = -0.01)) +
  NULL

legend <- get_legend(prop_plot + theme(legend.position = "right"))

prop_total <- plot_grid(total_plot ,
                        prop_plot ,
                        ncol = 1,
                        align = "vh")

prop_total_legend <- plot_grid(prop_total, ncol = 2, rel_widths = c(1, 0.4))

ggsave(prop_total, file = "plots/prop_total.jpg", dpi = 300, width = 7, height = 5, units = "in")


# Plot emergence ----------------------------------------------------------

(emerge_plot <- emerge_cond_posts %>% 
   ggplot(aes(x = date, color = trt, fill = trt)) +
   stat_pointinterval(position = position_dodge(width = 2), aes(y = .epred)) + 
   # geom_line(aes(group = interaction(trt, .draw)), alpha = 0.1) + 
   scale_color_colorblind() + 
   scale_fill_colorblind() +
   geom_point(data = emerge_reu_mg, aes(y = tot_mg_dm_m2_d, 
                                        fill = trt2), 
              position = position_jitterdodge(jitter.width = 0, dodge.width = 2), 
              alpha = 0.5, shape = 21, color = "black") +
   theme_tidybayes() +
   labs(subtitle = "a) Emerging insects (>99% chironomids)",
        y = expression('mg dry mass/'~m^2/d)) +
   theme(legend.title = element_blank(),
         axis.title.x =element_blank()) +
   geom_vline(xintercept=as.Date("2017-06-02"),linetype=2) +
   coord_cartesian(xlim = c(as.Date("2017-05-28"), as.Date("2017-06-29"))) +
   geom_segment(aes(x = as.Date("2017-06-05"), y = 1750, xend = as.Date("2017-06-02"), yend =1750),
                arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "black") +
   annotate(geom = "text", x = as.Date("2017-06-07"), y = 1750, size = 3.5, label = "Start of\nexperiment"))

# Plot benthics -----------------------------------------------------------

#total mg summed across all prey taxa
total_mg <- benthic_cond_posts %>% 
  group_by(.draw, date, trt) %>% 
  summarize(mg_dm = sum(.epred)) %>% 
  mutate(group = "All insects")

raw_mg <- ben_dm_tot %>% as_tibble() %>% 
  group_by(id, location, date, trt) %>% 
  summarize(mg_dm = sum(benthic_mg_dm_m2)) %>% 
  mutate(group = "All insects",
         date = mdy(date),
         trt = case_when(trt == "ctrl" ~ "fish", TRUE ~ "no fish"))

#total mg chiro only
total_chiro <- benthic_cond_posts %>% 
  filter(taxon == "chironomidae") %>% 
  group_by(.draw, date, trt) %>% 
  summarize(mg_dm = sum(.epred)) %>% 
  mutate(group = "Chironomidae")

raw_chiro <-  ben_dm_tot %>% as_tibble() %>% 
  filter(taxon == "chironomidae") %>% 
  group_by(id, location, date, trt) %>% 
  summarize(mg_dm = sum(benthic_mg_dm_m2)) %>% 
  mutate(group = "Chironomidae",
         date = mdy(date),
         trt = case_when(trt == "ctrl" ~ "fish", TRUE ~ "no fish"))

(benthic_plot <- total_mg %>% 
  ggplot(aes(x = date, y = mg_dm, color = trt, fill = trt)) +
  stat_pointinterval(position = position_dodge(width = 2)) + 
  # geom_line(aes(group = interaction(trt, .draw)), alpha = 0.1) + 
  scale_color_colorblind() + 
    scale_fill_colorblind() +
  geom_point(data = raw_mg, position = position_jitterdodge(jitter.width = 0, dodge.width = 2), 
             alpha = 0.5, shape = 21, color = "black") +
  theme_tidybayes() +
  labs(subtitle = "b) Benthic invertebrates (84-91% chironomids)",
       y = bquote('mg dry mass/'~m^2)) +
  theme(legend.title = element_blank(),
        axis.title.x =element_blank()) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2) +
  coord_cartesian(xlim = c(as.Date("2017-05-28"), as.Date("2017-06-29"))))

  




# Plot spiders -----------------------------------------------------------

spiders_cond_posts

(spider_plot <- spiders_cond_posts %>% 
    ggplot(aes(x = date, y = .epred, color = trt, fill = trt)) +
    stat_pointinterval(position = position_dodge(width = 2)) + 
    # geom_line(aes(group = interaction(trt, .draw)), alpha = 0.1) + 
    scale_color_colorblind() + 
    scale_fill_colorblind() +
    geom_point(data = spider_abund, aes(y = spid, fill = trt2), 
               position = position_jitterdodge(jitter.width = 0, dodge.width = 2), 
               alpha = 0.5, shape = 21, color = "black") +
    theme_tidybayes() +
    labs(subtitle = "c) Spider abundance",
         y = bquote('mg dry mass/'~m^2)) +
    theme(legend.title = element_blank(),
          axis.title.x =element_blank()) +
    geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
    coord_cartesian(xlim = c(as.Date("2017-05-28"), as.Date("2017-06-29"))))



# Combine plots -----------------------------------------------------------


emerge_ben_spid <- plot_grid(emerge_plot, benthic_plot, spider_plot, ncol = 1, align = "hv")
ggsave(emerge_ben_spid, file = "plots/emerge_ben_spid.tiff", dpi = 400, width = 6, height = 8)

