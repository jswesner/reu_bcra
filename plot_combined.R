library(tidyverse)
library(cowplot)
library(janitor)


# Load plots --------------------------------------------------------------

plot_benthic <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_benthic.rds?raw=true"))
plot_emerge <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_emerge.rds?raw=true"))
plot_spider <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_spider.rds?raw=true"))

emerge_plot_fishcage <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/emerge_plot_fishcage.rds?raw=true"))
benthic_plot_fishcage <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/benthic_plot_fishcage.rds?raw=true"))

plot_benthic_post_prior <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_benthic_post_prior.rds?raw=true"))
plot_emerge_post_prior <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_emerge_post_prior.rds?raw=true"))



# Combine plots - main figures -----------------------------------------------------------


a <- plot_emerge + theme(text = element_text(size = 14))+
ggtitle("a) Emerging insects (>99% chironomids)")
b <- plot_benthic + theme(text = element_text(size = 14))
c <- plot_spider + theme(text = element_text(size = 14))


plot_all <- plot_grid(a,b,c, ncol = 1,
          align = "v")

ggsave(plot_all, file = "plot_all.tiff", dpi = 600, width = 7, height = 9)
saveRDS(plot_all, file = "plot_all.rds")



# Combine plots - supplemental figures -----------------------------------------------------------
a2 <- benthic_plot_fishcage + theme(text = element_text(size = 14))
b2 <- emerge_plot_fishcage + theme(text = element_text(size = 14))

plot_all_fishcage <- plot_grid(a2,b2, ncol = 1,
                      align = "v")
ggsave(plot_all_fishcage, file = "plot_all_fishcage.tiff", dpi = 600, width = 7, height = 6)
saveRDS(plot_all_fishcage, file = "plot_all_fishcage.rds")




#prior v post
plot_prior_v_post <- plot_grid(plot_benthic_post_prior + guides(fill = F), 
                               plot_emerge_post_prior + theme(axis.title.y = element_blank()),
                               ncol = 2, align = "h")

ggsave(plot_prior_v_post, file = "plot_prior_v_post.tiff", dpi = 600, width = 12, height = 11 )





