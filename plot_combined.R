library(tidyverse)
library(cowplot)
library(janitor)


# Load plots --------------------------------------------------------------

plot_benthic <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_benthic.rds?raw=true"))
plot_emerge<- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_emerge.rds?raw=true"))
plot_spider<- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_spider.rds?raw=true"))

emerge_plot_fishcage <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/emerge_plot_fishcage.rds?raw=true"))
benthic_plot_fishcage <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/benthic_plot_fishcage.rds?raw=true"))

# Combine plots - main figures -----------------------------------------------------------

a <- plot_benthic + theme(text = element_text(size = 14))
b <- plot_emerge + theme(text = element_text(size = 14))
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
