rm(list = ls())

# load libraries ----------------------------------------------------------
library(tidyverse)
library(xlsx)

# read R function ---------------------------------------------------------
wdir <- "./chapters/chapter4/"
source(paste0(wdir, "scripts/func/read_DDD_results.R"))

# set parameters ----------------------------------------------------------
root.dir = "U:/Lettuce_RTD_analysis" # root files directory
df <- read_csv(paste(root.dir, "3D_output_dirs.csv", sep = "/"))
dataset <- unique(df$dataset)
ref.idx <- unique(df$ref.idx)

sig.DE.genes <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "Significant DE genes list and statistics.csv"
)

sig.DAS.genes <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "Significant DAS genes list and statistics.csv"
)

df <- left_join(
  x = sig.DAS.genes,
  y = sig.DE.genes,
  by = join_by(target == target,
               contrast == contrast),
  keep = TRUE
)

df <- df |> 
  mutate(DD = ifelse(is.na(up.down.y),  "DASonly", "DE&DAS")) |>
  mutate(up.down.y = ifelse(is.na(up.down.y), "DASonly", up.down.y)) |> 
  mutate(hpi = as.numeric(str_extract(contrast.x, "\\d+\\.?\\d*"))) |>
  group_by(hpi, up.down.y) |>
  summarise(count = n())

df$up.down.y <- factor(df$up.down.y, levels = c("down-regulated", "up-regulated", "DASonly")) 

df <- df |> 
  ungroup() |>
  group_by(hpi) |>
  arrange(hpi, up.down.y) |> 
  mutate(label_y = cumsum(count))
  
df$up.down.y <- factor(df$up.down.y, levels = c("DASonly", "up-regulated", "down-regulated")) 

plt <- ggplot(data = df, mapping =  aes(x = hpi, y = count, fill = up.down.y, color = up.down.y)) +
  geom_col(position = "stack") +
  geom_text(aes(y = label_y, label = count),
            vjust = -0.5, color = "black") +
  scale_fill_manual(
    values = c(
      "up-regulated" = "#E69F00",
      "down-regulated" = "#56B4E9",
      "DASonly" = "darkgrey"
    ),
    labels = c(
      "up-regulated" = "DAS & DE\nup-regulated",
      "down-regulated" = "DAS & DE\ndown-regulated",
      "DASonly" = "only DAS"
    )
  ) +
  scale_color_manual(
    values = c(
      "up-regulated" = "#E69F00",
      "down-regulated" = "#56B4E9",
      "DASonly" = "darkgrey"
    ),
    labels = c(
      "up-regulated" = "DAS & DE\nup-regulated",
      "down-regulated" = "DAS & DE\ndown-regulated",
      "DASonly" = "only DAS"
    )
  ) +
  scale_x_continuous(breaks = seq(9, 48, 3),
                     expand = c(0, 1)) +
  scale_y_continuous(
    breaks = seq(0, 1800, 200),
    limits = c(0, 1800),
    expand = c(0, 1)
  ) +
  labs(x = 'Hours post inoculation',
       y = 'Number of DAS genes') +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

plt
# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_DAS_genes.png"),
  width = 6,
  height = 4.5,
  units = "in",
  res = 300
)
plt
dev.off()

svg(
  file = paste0(wdir, "figures/time-series_DAS_genes.svg"),
  width = 6,
  height = 4.5,
)
plt
dev.off()

# DTU transcripts
# set parameters ----------------------------------------------------------
root.dir = "U:/Lettuce_RTD_analysis" # root files directory
df <- read_csv(paste(root.dir, "3D_output_dirs.csv", sep = "/"))
dataset <- unique(df$dataset)
ref.idx <- unique(df$ref.idx)
sig.DE.tx <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "Significant DE transcripts list and statistics.csv"
)

sig.DTU.tx <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "Significant DTU transcripts list and statistics.csv"
)

df <- left_join(
  x = sig.DTU.tx,
  y = sig.DE.tx,
  by = join_by(target == target,
               contrast == contrast),
  keep = TRUE
)

df <- df |> 
  mutate(DD = ifelse(is.na(up.down.y),  "DTUonly", "DE&DTU")) |>
  mutate(up.down.y = ifelse(is.na(up.down.y), "DTUonly", up.down.y)) |> 
  mutate(hpi = as.numeric(str_extract(contrast.x, "\\d+\\.?\\d*"))) |>
  group_by(hpi, up.down.y) |>
  summarise(count = n())

df$up.down.y <- factor(df$up.down.y, levels = c("down-regulated", "up-regulated", "DTUonly")) 
levels(df$up.down.y)
df <- df |> 
  ungroup() |>
  group_by(hpi) |>
  arrange(hpi, up.down.y) |> 
  mutate(label_y = cumsum(count))

df$up.down.y <- factor(df$up.down.y, levels = c("DTUonly", "up-regulated", "down-regulated")) 

plt <- ggplot(data = df, mapping =  aes(x = hpi, y = count, fill = up.down.y, color = up.down.y)) +
  geom_col(position = "stack") +
  geom_text(aes(y = label_y, label = count),
            vjust = -0.5, color = "black") +
  scale_fill_manual(
    values = c(
      "up-regulated" = "#E69F00",
      "down-regulated" = "#56B4E9",
      "DTUonly" = "darkgrey"
    ),
    labels = c(
      "up-regulated" = "DTU & DE\nup-regulated",
      "down-regulated" = "DTU & DE\ndown-regulated",
      "DTUonly" = "only DTU"
    )
  ) +
  scale_color_manual(
    values = c(
      "up-regulated" = "#E69F00",
      "down-regulated" = "#56B4E9",
      "DTUonly" = "darkgrey"
    ),
    labels = c(
      "up-regulated" = "DTU & DE\nup-regulated",
      "down-regulated" = "DTU & DE\ndown-regulated",
      "DTUonly" = "only DTU"
    )
  ) +
  scale_x_continuous(breaks = seq(9, 48, 3),
                     expand = c(0, 1)) +
  scale_y_continuous(
    breaks = seq(0, 2200, 200),
    limits = c(0, 2200),
    expand = c(0, 1)
  ) +
  labs(x = 'Hours post inoculation',
       y = 'Number of DTU transcripts') +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_DTU_tx.png"),
  width = 6,
  height = 4.5,
  units = "in",
  res = 300
)
plt
dev.off()

svg(
  file = paste0(wdir, "figures/time-series_DTU_tx.svg"),
  width = 6,
  height = 4.5,
)
plt
dev.off()
