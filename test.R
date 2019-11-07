diet_mgdm %>% 
  filter(species == "quillback")
  filter(prey_taxon == "fish_yoy",
         mg_diet_dm>0) %>% 
  print(n = Inf)

post1000_agg_0 <- post1000_agg %>% 
  left_join(d_0) %>% 
  ungroup()


post1000_agg_0 %>%
  filter(prey_taxon != "fish_yoy") %>% 
  mutate(mg_dm_diet_corrected = mg_dm_diet*correction) %>% 
  ggplot(aes(x = reorder(prey_taxon,mg_dm_diet), y = mg_dm_diet_corrected, 
             fill = date2,
             group = interaction(prey_taxon, date2))) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(-.5))+
  ylim(c(0,30)) +
  facet_grid(.~species) +
  #geom_point(alpha = 0.5) +
  #geom_violin() + +
  #scale_y_log10()  +
  scale_fill_grey(name = "sample date") +
  #scale_color_grey() +
  coord_flip() +
  theme_classic() +
  xlab("Prey taxon") +
  ylab("Diet (mg dry mass)") 


diet_mgmd %>% 
  group_by(prey_taxon) %>% 
  summarize(max = max(mg_diet_dm)) %>% 
  arrange(max)
  ggplot(aes(x = date, y = mg_diet_dm + 0.01, color = prey_taxon)) +
  geom_point(position = position_jitter(), alpha =0.5) +
  facet_wrap(~species, scales = "free")+
  scale_y_log10()

  
d_0 %>% 
  group_by(prey_taxon) %>% 
  summarize(max = max(correction))

  