rm(list = ls())
# Gene expression profile of DE novel genes
# load libraries ----------------------------------------------------------
require(tidyverse)
# genes of interest -------------------------------------------------------
gene <- c(
  "LSATv11_C01_002607",
  "LSATv11_C04_023408",
  "LSATv11_C06_045243",
  "LSATv11_C03_020086",
  "LSATv11_C04_023358",
  "LSATv11_C05_038176"
)

# source plotting function ------------------------------------------------
wdir <- "./chapters/chapter4/"
source(paste0(wdir, "scripts/func/plot_timeseries.R"))
# plot expression profile of each gene ------------------------------------
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
                     labels = NULL,
                     nrow = 3,
                     ncol = 1)



right.panel <-
  cowplot::plot_grid(plt4,
                     plt5,
                     plt6,
                     labels = NULL,
                     nrow = 3,
                     ncol = 1)

plt <- cowplot::plot_grid(
  left.panel,
  right.panel,
  labels = "AUTO",
  nrow = 1,
  ncol = 2
)

# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_novel_DE_genes_exp.png"),
  width = 9,
  height = 12,
  units = "in",
  res = 300
)
plt
dev.off()

svg(
  file = paste0(wdir, "figures/time-series_novel_DE_genes_exp.svg"),
  width = 9,
  height = 12,
)
plt
dev.off()
