rm(list = ls())
# Differentially expressed lncRNAs in the time-series data
require(tidyverse)
## read data --------------------------------------------------------------
wdir <- "./chapters/chapter4/"

df <- readxl::read_xlsx(
  path = paste0(wdir, "suppl/suppl_2.xlsx"),
  sheet = "Table S6",
  skip = 4
)
## tidy-up dataframe -------------------------------------------------------
df <- df |>
  mutate(Direction = factor(Direction, levels = c("up-regulated", "down-regulated"))) |>
  group_by(HPI) |>
  mutate(label_y = cumsum(N_DE_Genes))

y_max = 100 * (max(df$label_y) %/% 100 + 1) # maximum limit of y axis
## plot DE lncRNA genes ----------------------------------------------------
plt <-
  ggplot(data = df, aes(x = HPI, y = N_DE_Genes, fill = Direction)) +
  geom_col(position = "stack") +
  geom_text(
    aes(y = label_y, label = N_DE_Genes),
    colour = "black",
    hjust = 1.1,
    angle = 90
  ) +
  scale_fill_manual(values = c(
    "up-regulated" = "#E69F00",
    "down-regulated" = "#56B4E9"
  )) +
  scale_x_continuous(
    name = "Hours post inoculation",
    limits = c(9, 48),
    breaks = seq(9, 48, by = 3),
    oob = scales::squish_infinite
  ) +
  scale_y_continuous(
    name = "Number of DE lncRNA genes",
    limits = c(0, y_max),
    breaks = seq(0, y_max, by = 100),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank()
  )
plt
# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_DE_lncRNA_genes.png"),
  width = 6,
  height = 4.5,
  units = "in",
  res = 300
)
plt
dev.off()

svg(
  file = paste0(wdir, "figures/time-series_DE_lncRNA_genes.svg"),
  width = 6,
  height = 4.5,
)
plt
dev.off()

