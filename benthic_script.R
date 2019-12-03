library(tidyverse)
library(lubridate)
library(brms)
library(ggridges)
library(RCurl)
library(janitor)



# Load data ---------------------------------------------------------------
#benthic abundance (insects only)
ben <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/benthic62617.csv")) %>% 
  clean_names() %>% 
  mutate(chironomidae = chiro_larvae + pupae) %>% 
  select(-chiro_larvae, -pupae)

#benthic dry mass of individuals
ben_dm_ind <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/benthic_drymass_141517.csv")) %>% 
  filter(mg_dm_per_ind > 0) %>% 
  mutate(taxon = str_replace(taxon, "Baetidae", "ephemeroptera"),
         taxon = str_replace(taxon, "Caenidae", "ephemeroptera"))

#Bayesian brms model of length-weight regression
brm_ind_mg <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/benthic_ind_mg.RDS?raw=true"))

#Bayesian brms model of benthic dry mass
brm_ben_m2 <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/brm_ben_m2.RDS?raw=true"))


#Bayesian model of benthic individual dry mass
#brm_ind_mg <- brm(mg_dm_per_ind ~ 0 + taxon + (1|experiment), data = ben_dm_ind,
#                  family = Gamma(link = "log"),
#                  prior = c(prior(normal(0,1), class = "b")),
#                 cores = 4)

#check model
# brm_ind_mg
# pp_check(brm_ind_mg, type = "boxplot")

#saveRDS(brm_ind_mg, file = "benthic_ind_mg.RDS")

#extract posterior and estimate mean
ind_mg_post <- posterior_samples(brm_ind_mg) %>% 
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


ggsave(plot_ben_cage_v_amb, file = "plot_ben_cage_v_amb.tiff", dpi = 600, width = 6.5, height =4.5)



#make data to analyze - total benthic mg
ben_dm_tot <- merge(ben_clean,ind_mg_mean) %>% 
  mutate(sample_m2 = 0.76*0.3,
         benthic_mg_dm_m2 = (abundance*mean_mgdm)/sample_m2,
         benthic_mg01 = 0.01 + benthic_mg_dm_m2) 




# Bayesian model ----------------------------------------------------------


#brm_ben_m2 <- brm(benthic_mg01~date*trt*taxon + (1|location), data=ben_dm_tot,family=Gamma(link="log"),
#prior=c(prior(normal(0,2),class="Intercept"), 
#prior(normal(0,1),class="b"),
#prior(cauchy(0,1),class = "sd")))

#check model
brm_ben_m2
pp_check(brm_ben_m2, type = "boxplot")

#saveRDS(brm_ben_m2, file = "brm_ben_m2.RDS")

# Extract conditional posteriors ------------------------------------------

date = unique(ben_dm_tot$date)
trt = unique(ben_dm_tot$trt)
taxon = unique(ben_dm_tot$taxon)
newdata = expand.grid(date,trt,taxon) %>% 
  rename(date = Var1,
         trt = Var2,
         taxon = Var3)

marg_post <- fitted(brm_ben_m2, summary = F, newdata = newdata, re_formula = NA)
names = newdata %>% unite("colnames", date:taxon)
colnames(marg_post) <- names$colnames

marg_post2 <- marg_post %>% as_tibble(marg_post)%>%
  mutate(iter = 1:nrow(marg_post)) %>% 
  gather(key,mg_dm_m2, -iter)%>%
  separate(key, c("date","trt","taxon"), sep = "_") %>% 
  mutate(date = mdy(as.factor(date)))

total_mg <- marg_post2 %>% 
  group_by(iter, date, trt) %>% 
  summarize(mg_dm = sum(mg_dm_m2))


# Plot posterior ----------------------------------------------------------
raw_ben_plot <- as_tibble(ben_dm_tot) %>% 
  mutate(date = mdy(as.factor(date)),
         trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, id, location, trt) %>% 
  summarize(mg_dm = sum(benthic_mg_dm_m2))


plot_benthic <- total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                       "exc" = "no fish"))) %>%
  ggplot(aes(x = date, y = mg_dm, fill = trt)) +
  geom_boxplot(aes(group = interaction(trt, date)), outlier.shape = NA,
               position = position_dodge(width = 2),
               width = 1.5)+
  theme_classic()+
  scale_fill_grey(start = 0.9, end = 0.4)+
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  scale_y_continuous(limits = c(0, 4800),
                     breaks = c(0,1000, 2000, 3000, 4000, 5000)) +
  theme(legend.title = element_blank(),
        axis.title.x =element_blank(),
        text = element_text(size = 20))+
  geom_point(data = raw_ben_plot, aes(y = mg_dm, x = date, fill = trt), 
             position = position_dodge(width = 2),
             shape = 16, size = 1.3) +
  ylab(bquote('mg dry mass/'~m^2)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  #scale_y_log10()+
  ggtitle("b) Benthic insects")


ggsave(plot_benthic, file = "plot_benthic.tiff", dpi = 600, width = 7, height = 3.5, units = "in")
saveRDS(plot_benthic, file = "plot_benthic.rds")


# Summarize posterior -----------------------------------------------------
#summary stats by date and treatment
tot_benthic <- total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no_fish"))) %>% 
  group_by(date, trt) %>% 
  summarize(mean = mean(mg_dm),
            median = median(mg_dm),
            sd = sd(mg_dm),
            low95 = quantile(mg_dm, probs = 0.025),
            high95 = quantile(mg_dm, probs = 0.975))


write.csv(tot_benthic, file = "tot_benthic.csv")


#proportion difference by trt
total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, trt) %>% 
  pivot_wider(names_from = trt,
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(prop_reduction = 1- fish/no_fish) %>% 
  summarize(mean = mean(prop_reduction),
            median = median(prop_reduction),
            sd = sd(prop_reduction),
            low95 = quantile(prop_reduction, probs = 0.025),
            high95 = quantile(prop_reduction, probs = 0.975))

#probability of a difference
total_mg %>%
  ungroup() %>% 
  mutate(trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish"))) %>% 
  group_by(date, trt) %>% 
  pivot_wider(names_from = trt,
              values_from = mg_dm) %>% 
  clean_names() %>% 
  mutate(diff = fish-no_fish) %>% 
  summarize(prob_nofish_higher = sum(diff<0)/4000)




# Fish in cages -----------------------------------------------------------

fish_in_cages <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/fish_in_cages.csv")) %>% 
  mutate(date = mdy(date))

fish_in_cages_select <- fish_in_cages %>%  select(fish_in_cages, station)

benthic_plot_fishcage <- raw_ben_plot %>% 
  unite("station", c("id","location")) %>%
  left_join(fish_in_cages_select) %>% 
  mutate(fish_in_cages = ifelse(date == "2017-06-02",0, fish_in_cages)) %>% 
  ggplot(aes(x = date, y = mg_dm, fill = trt)) +
  geom_point(position = position_dodge(width = 2),
             shape = 21, alpha = 0.8, aes(size = fish_in_cages)) +
  scale_fill_grey(name = "")+
  scale_size_continuous(name = "# fish in cages on 2017-06-30")+
  theme_classic()+
  theme(axis.title.x =element_blank())+
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  ylab(bquote('mg dry mass/'~m^2)) +
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  annotate("text",x=as.Date("2017-06-02")+7.5,y=3200,label="start of experiment")+
  geom_segment(aes(x = as.Date("2017-06-04"), y = 3200, 
                   xend=as.Date("2017-06-02"), yend = 3200),
               size = 0.1)+
  #scale_y_log10()+
  ggtitle("a) Benthic insects")+
  #annotate("text",x=as.Date("2017-06-02")+7.5,y=320,label="start of experiment")+
  #geom_segment(aes(x = as.Date("2017-06-05"), y = 320, xend=as.Date("2017-06-02"), yend = 320))+
  #scale_y_log10()+
  NULL



ggsave(benthic_plot_fishcage, file = "benthic_plot_fishcage.tiff", dpi = 600, width = 7, height = 3.5, units = "in")
saveRDS(benthic_plot_fishcage, file = "benthic_plot_fishcage.rds")

