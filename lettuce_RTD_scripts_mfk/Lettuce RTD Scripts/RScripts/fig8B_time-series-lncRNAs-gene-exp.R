rm(list = ls())
# Gene expression profile of DE lncRNAs
# load libraries ----------------------------------------------------------
library(tidyverse)
# genes of interest -------------------------------------------------------
gene <- c(
  "LSATv11_C03_017500",
  "LSATv11_C05_036677",
  "LSATv11_C09_061369",
  "LSATv11_C02_013893",
  "LSATv11_C04_022579",
  "LSATv11_C07_048246"
)

# source plotting function ------------------------------------------------
wdir <- "./chapters/chapter4/"
source(paste0(wdir, "scripts/func/plot_lncRNA_timeseries.R"))
## plot expression profile of new genes ------------------------------------
plt1 <- plot.timeseries(paste0(gene[1]))
plt2 <- plot.timeseries(paste0(gene[2]))
plt3 <- plot.timeseries(paste0(gene[3]))
plt4 <- plot.timeseries(paste0(gene[4]))
plt5 <- plot.timeseries(paste0(gene[5]))
plt6 <- plot.timeseries(paste0(gene[6]))

# save multiple plots into a grid -----------------------------------------
left.panel <-
  cowplot::plot_grid(plt1,
                     plt2,
                     plt3,
                     nrow = 3,
                     ncol = 1)



right.panel <-
  cowplot::plot_grid(plt4,
                     plt5,
                     plt6,
                     nrow = 3,
                     ncol = 1)

plt <- cowplot::plot_grid(
  left.panel,
  right.panel,
  nrow = 1,
  ncol = 2
)

plt

# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_lncRNAs_gene_exp.png"),
  width = 9,
  height = 12,
  units = "in",
  res = 300
)
plt
dev.off()

svg(
  file = paste0(wdir, "figures/time-series_lncRNAs_gene_exp.svg"),
  width = 9,
  height = 12,
)
plt
dev.off()
