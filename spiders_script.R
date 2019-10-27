#Packages
library(tidyverse)
library(lubridate)
library(brms)
library(ggridges)
library(RCurl)
library(janitor)

#Summary of counts within traps

# Load data ---------------------------------------------------------------

#spider abundance
spider_abund <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/spiders_dens.csv"))
emerge_reu_mg <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/reu_bcra/master/emerge_reu_mg.csv"))

#Bayesian brms spider abundance model
m2 <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/spider_brms.RDS?raw=true"))


#emergence data for dates near spider collections
emerge_reu_mg_lasttwodates <- as_tibble(emerge_reu_mg) %>% 
  filter(date2=="2017-06-16"|date2=="2017-06-23") %>% 
  mutate(spider_emerge_date = ifelse(date2=="2017-06-16",1,2),
         date3 = as.Date(date2),
         spider_emerge_date = as.Date(date3))%>% 
  select(date3,stat,stat2,loc,trt,tot_mg_dm_m2_d,spider_emerge_date) %>% 
  filter(trt!="Ambient")

#spider data to merge with emergence
spider_abund_tomerge <- as_tibble(spider_abund) %>% 
  mutate(spider_emerge_date = ymd(as.Date(date,"%m/%d/%Y")),
         date2 = ymd(as.Date(date,"%m/%d/%Y"))) %>% 
  select(date2,stat,stat2,loc,trt,spid,spider_emerge_date) %>% 
  filter(trt!="Ambient")

spider_emerge <- spider_abund_tomerge %>% 
  left_join(emerge_reu_mg_lasttwodates,by="stat2")


spider_abund2 <- as_tibble(spider_abund) %>% 
  filter(trt!="Ambient") %>% 
  mutate(date = as.factor(ymd(as.Date(date,"%m/%d/%Y"))),
         trt = ifelse(trt=="Cage Control","ctrl",
                      ifelse(trt=="exc","exc","exc"))) 



#Bayesian model
spider_abund2$spid01 <- 0.01+spider_abund2$spid
#m2<-brm(spid~trt*date,data=spider_abund2,family=poisson(link="log"),
#        prior=c(prior(normal(0,100),class="b"),
#                prior(normal(0,20),class="Intercept")),
#        iter = 3000, chains = 4, cores=4)

#check model
m2
pp_check(m2, type = "boxplot")

#saveRDS(m2, file = "spider_brms.RDS")

# Extract conditional posteriors ------------------------------------------

marg_spiders <-marginal_effects(m2,method="fitted",effects="date:trt")
marg_spiders_fit <- fitted(m2, newdata=marg_spiders$`date:trt`,summary=F)
columns_spiders<- paste(marg_spiders$`date:trt`$date,"_",marg_spiders$`date:trt`$trt)
colnames(marg_spiders_fit) <- columns_spiders
marg_spiders_fit <- as.data.frame(marg_spiders_fit)
marg_spiders_fit$iter <- 1:nrow(marg_spiders_fit)

marg_spiders_fit_plot <- as_tibble(marg_spiders_fit) %>% 
  gather(date_trt, spider_abundance,-iter) %>%
  separate(date_trt, c("date","trt"), sep=' _ ')


# Plot posteriors ---------------------------------------------------------
raw_spider_plot <- spider_abund2 %>% 
  mutate(date = ymd(date),
         trt = str_replace_all(trt,c("ctrl" = "fish",
                                     "exc" = "no fish")),
         spider_abundance = spid)
  

plot_spider <- marg_spiders_fit_plot %>% 
  mutate(date = ymd(date),
         trt = str_replace_all(trt,c("ctrl" = "fish",
                                       "exc" = "no fish"))) %>%
  ggplot(aes(x = date, y = spider_abundance, fill = trt)) +
  geom_boxplot(aes(group = interaction(trt, date)), outlier.shape = NA,
               position = position_dodge(width = 2),
               width = 1.5)+
  theme_classic()+
  scale_fill_grey(start = 0.9, end = 0.4)+
  theme(legend.title = element_blank(),
        axis.title.x =element_blank(),
        text = element_text(size = 20))+
  coord_cartesian(xlim = as.Date(c('2017-05-28', '2017-06-29'), 
                                 format="%Y-%m-%d")) +
  geom_point(data = raw_spider_plot, aes(fill = trt), 
             position = position_dodge(width = 2),
             shape = 21, size = 1.3) +
  ylab(expression(paste("Spider abundance (#/cage)")))+
  geom_vline(xintercept=as.Date("2017-06-02"),linetype=2)+
  #annotate("text",x=as.Date("2017-06-02")+7.5,y=320,label="start of experiment")+
  #geom_segment(aes(x = as.Date("2017-06-05"), y = 320, xend=as.Date("2017-06-02"), yend = 320))+
  #scale_y_log10()+
  ggtitle("c) Spider abundance")

ggsave(plot_spider, file = "plot_spider.tiff", dpi = 600, width = 7, height = 3.5, units = "in")

  

# Summarize posterior -----------------------------------------------------

#mean spiders
marg_spiders_fit_plot %>% 
  group_by(date,trt) %>% 
  summarize(mean = mean(spider_abundance),
            median = median(spider_abundance),
            sd = sd(spider_abundance),
            low95 = quantile(spider_abundance,probs=0.025),
            high95 = quantile(spider_abundance,probs=0.975)) %>% 
  mutate_if(is.numeric,round,1) 

#differences in treatments
marg_spiders_fit_plot %>% 
  group_by(date,trt) %>% 
  spread(trt,spider_abundance) %>% 
  mutate(diff = exc-ctrl) %>% 
  summarize(low95 = quantile(diff,probs=0.025),
            median = median(diff),
            high95 = quantile(diff,probs=0.975)) %>% 
  mutate_if(is.numeric,round,1)

#probability of difference
marg_spiders_fit_plot %>% 
  spread(trt,spider_abundance) %>% 
  group_by(date) %>% 
  mutate(diff = exc-ctrl) %>% 
  summarize(prob_exc_more_than_ctrl = sum(diff>0)/4000) %>% 
  mutate_if(is.numeric,round,2)

#Regression across all tanks


ggplot(spider_emerge,aes(x=tot_mg_dm_m2_d,y=spid))+
  geom_point(alpha=0.5, 
             position = position_jitter(width = 0, height = 0.02))+
  geom_smooth(method="lm")+
  facet_wrap(~spider_emerge_date.x)

