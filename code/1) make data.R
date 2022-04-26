#load packages
source("src/packages.R")

# Diet --------------------------------------------------------------------
# Load data  
diet_mgdm <- read_csv("data/raw_data/diet_mgdm.csv")

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
                                                    "mayfly" = "ephemeroptera"))) %>% 
  separate(prey_taxon, c("prey_taxon", "prey_stage")) %>% 
  ungroup() %>% 
  mutate(prey_stage = str_trim(case_when(prey_stage == "ad" ~ "a",
                                         is.na(prey_stage) ~ "unknown",
                                         TRUE ~ prey_stage))) %>% 
  unite("prey_taxon", prey_taxon:prey_stage, sep = "_") %>% 
  mutate(prey_mg_inverse = 1/prey_mg_dm)

write_csv(d, file = "data/derived_data/d.csv")

d_tot <- d %>% 
  separate(prey_taxon, c("prey_taxon", "prey_stage")) %>% 
  mutate(prey_stage = case_when(prey_stage %in% c("a", "p") ~ "pa", TRUE ~ "not_pa")) %>% 
  group_by(id, species, date2, prey_stage) %>% 
  summarize(mg_dmtotal = sum(mg_diet_dm)) %>% 
  arrange(-mg_dmtotal)

write_csv(d_tot, file = "data/derived_data/d_tot.csv")

d_groups <- d %>% 
  group_by(prey_taxon) %>% 
  summarize(total = sum(mg_diet_dm)) %>% 
  mutate(global = sum(total),
         proportion = total/global) %>% 
  arrange(-proportion) %>% 
  mutate(prey_taxon_grouped = case_when(proportion >= 0.001 | prey_taxon == "coleo_a" ~ prey_taxon,
                                        TRUE ~ "other_unknown"))

d_grouped <- d %>% left_join(d_groups %>% select(prey_taxon, prey_taxon_grouped)) %>%
  group_by(date2, species, prey_taxon_grouped, id) %>% 
  summarize(mg_diet_dm = sum(mg_diet_dm),
            mg_diet_dm01 = mg_diet_dm + 0.1)

rm(d_groups)

# fish community
fish_community <- read_csv("data/raw_data/fish_community.csv") %>% 
  mutate(date2 = case_when(date == as.Date("2017-06-12") ~ as.Date("2017-06-13"), 
                           TRUE ~ date),
         sampling = "density")


fish_totals <- fish_community %>% 
  group_by(date2, species) %>% 
  summarize(total_abund_raw = sum(abund)) %>% 
  group_by(species) %>% 
  mutate(average = mean(total_abund_raw),
         total_abund = case_when(total_abund_raw == 0 ~ average, TRUE ~ total_abund_raw)) %>% 
  filter(date2 != as.Date("2017-05-31")) %>% 
  mutate(date2 = as.factor(date2),
         species = case_when(species == "river_shiner" ~ "rivershiner",
                             species == "spotfin_shiner" ~ "spotfin",
                             TRUE ~ species)) 

write_csv(fish_totals, file = "data/raw_data/fish_totals.csv")

# Emergence ---------------------------------------------------------------
# load data
emerge_reu <- read_csv("data/raw_data/emerge_reu.csv")
chiro_ind <- read_csv("data/raw_data/chiro_ind.csv")  #chiro_individual_lengths

# Length-mass regression - chiros

#plot to check for differences in chiro length over time and treatments
chiro_ind %>%
  mutate(date = mdy(date)) %>%
  ggplot(aes(x=date,y=mm,color=trt,group=interaction(trt,date)))+
  geom_violin(position=position_dodge(width=8))+
  geom_point(position=position_dodge(width=8),size=3,alpha=0.2)+
  #scale_y_log10()+
  theme_classic()+
  scale_color_grey()+
  ggtitle("Violin plot of chironomid lengths (mm)")+
  NULL

#add parameters (a,b) to convert length to mass
chiro_mg <- chiro_ind %>%
  mutate(a = 0.1, #parameter
         b = 1.57, #parameter
         param_source = "Sabo et al. 2002 (Table 2)", #source
         param_taxon = "Nematocera") %>% #source taxon
  mutate(mg_dm = a*mm^b) #generate biomass

write_csv(chiro_mg, file = "data/derived_data/chiro_mg.csv")

# Fit model of chironomid average dry mass 
chiro_mg_dm <- brm(mg_dm ~ 1, data=chiro_mg, family=Gamma(link="log"),
                prior=prior(normal(0,2),class="Intercept"),
                file = "models/emerge_chiro_mg_dm.rds", 
                file_refit = "on_change")

#check model
chiro_mg_dm
pp_check(chiro_mg_dm, type = "boxplot")

#get mean and sd of chiro mass to sample from
as_draws_df(chiro_mg_dm) %>% 
  summarize(mean = mean(exp(b_Intercept)))

#Convert emergence abundance to dry mass (mg)

#Chironomid dry mass based on model above using direct measures of mm
#Dry mass of other taxa were estimated from previous collections at the same site or nearby sites (Wesner et al. in review)

set.seed(2020)
emerge_reu_mg <- emerge_reu %>%
  replace(is.na(.), 0) %>%
  mutate(chiro_mg_dm = chir*.99,
         doli_mg_dm = doli*.65,
         tric_mg_dm = tric*.8,
         odo_mg_dm = odo*6.2,
         ephem_mg_dm = ephe*1.1,
         tot_mg_dm = chiro_mg_dm+doli_mg_dm+tric_mg_dm+odo_mg_dm+ephem_mg_dm,
         tot_mg_dm_m2_d = tot_mg_dm/area/days,
         id = paste0(stat,"_",loc))

write_csv(emerge_reu_mg, file = "data/derived_data/emerge_reu_mg.csv")  


# Benthics ----------------------------------------------------------------
# the code below converts benthic abundance and length data to benthic dry mass
#benthic abundance (insects only)
ben <- read_csv("data/raw_data/benthic62617.csv") %>% 
  clean_names() %>% 
  mutate(chironomidae = chiro_larvae + pupae) %>% 
  select(-chiro_larvae, -pupae)

#benthic dry mass of individuals
ben_dm_ind <- read_csv("data/raw_data/benthic_drymass141517.csv") %>% 
  filter(mg_dm_per_ind > 0) %>% 
  mutate(taxon = str_replace(taxon, "Baetidae", "ephemeroptera"),
         taxon = str_replace(taxon, "Caenidae", "ephemeroptera"))

#Bayesian model of benthic individual dry mass
brm_ind_mg <- brm(mg_dm_per_ind ~ 0 + taxon + (1|experiment), 
                  data = ben_dm_ind,
                 family = Gamma(link = "log"),
                 prior = c(prior(normal(0,1), class = "b")),
                cores = 4,
                file = "models/brm_ind_mg.RDS",
                file_refit = "on_change")

#check model
brm_ind_mg
pp_check(brm_ind_mg, type = "boxplot")


#extract posterior and estimate mean
ind_mg_post <- as_draws_df(brm_ind_mg) %>% 
  select(starts_with('b')) %>% 
  gather()

ind_mg_mean <- ind_mg_post %>% 
  group_by(key) %>% 
  summarize(mean_mgdm = mean(exp(value))) %>% 
  mutate(taxon = tolower(str_remove(key, "b_taxon"))) %>% 
  select(-key) %>% 
  add_row(taxon="tipulidae",mean_mgdm=1.4) %>% 
  mutate(taxon = ifelse(taxon == "worms", "olgigochaetes", taxon))
#replace Caenidae with Ephem. #most mayflies are Caenids anyway at this 
#site, so this allows easer matching of names later on.
#Tipulidae were not measured for length. Instead, the mean of tipulid 
#lengths and regression to mass from Benke et al. (1999) was used to 
#estimate Tipulid length. Thus, there is no variation around this estimate
#(or at least it's very small - but non-zero to allow the gamma parameters
#to be fit).


#Convert abundance to biomass. First, clean benthic data to include only benthic organisms that emerged, and remove previous summaries. Then gather by taxon and add the dry mass data

ben_clean <- ben %>%
  select(-notes,-total,-daphnia,-isopodia,-hemiptera,-copepod,-snails,-oligochaetes) %>%
  gather(taxon,abundance,"ceratopogonidae":"chironomidae")

#plot ambient vs cage
plot_ben_cage_v_amb <- ben_clean %>% 
  mutate(date = mdy(date)) %>% 
  filter(trt2 != "Exclusion",
         date != "2017-06-02") %>% 
  ggplot(aes(x = reorder(taxon, -abundance), y = abundance, color = trt2, shape = trt2))+
  geom_point(size = 2,position = position_jitterdodge(dodge.width = 0.6,
                                                      jitter.width = 0),
             alpha = .8)+
  theme_bw()+
  facet_grid(date~., scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 10))+
  scale_y_log10() +
  scale_color_grey(start = 0.2, end = 0.7)+
  theme(legend.title = element_blank())+
  xlab("prey_taxon")+
  ylab("Number per benthic sample")+
  #coord_flip() +
  NULL


ggsave(plot_ben_cage_v_amb, file = "plots/plot_ben_cage_v_amb.tiff", dpi = 600, width = 6.5, height =4.5)

#make data to analyze - total benthic mg
ben_dm_tot <- merge(ben_clean,ind_mg_mean) %>% 
  mutate(sample_m2 = 0.76*0.3,
         benthic_mg_dm_m2 = (abundance*mean_mgdm)/sample_m2,
         benthic_mg01 = 0.01 + benthic_mg_dm_m2,
         benthic_g01 = benthic_mg_dm_m2/1000 + 0.001)

write_csv(ben_dm_tot, file = "data/derived_data/ben_dm_tot.csv")

rm(ben)
rm(ben_dm_ind)
rm(ben_clean)
rm(brm_ind_mg)
rm(d_grouped)
rm(d)
rm(diet_mgdm)
rm(ind_mg_mean)
rm(ind_mg_post)
rm(plot_ben_cage_v_amb)
