library(tidyverse)
library(cowplot)
library(janitor)


# Load plots --------------------------------------------------------------

plot_benthic <- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_benthic.rds?raw=true"))
plot_emerge<- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_emerge.rds?raw=true"))
plot_spider<- readRDS(url("https://github.com/jswesner/reu_bcra/blob/master/plot_spider.rds?raw=true"))


# Combine plots -----------------------------------------------------------

a <- plot_benthic + theme(text = element_text(size = 14))
b <- plot_emerge + theme(text = element_text(size = 14))
c <- plot_spider + theme(text = element_text(size = 14))


plot_all <- plot_grid(a,b,c, ncol = 1,
          align = "v")

ggsave(plot_all, file = "plot_all.tiff", dpi = 600, width = 7, height = 9)
saveRDS(plot_all, file = "plot_all.rds")
