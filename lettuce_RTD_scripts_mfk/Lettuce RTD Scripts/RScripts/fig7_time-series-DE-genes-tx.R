rm(list = ls())
# DE gene and transcripts over the course of B. cinerea infection
require(tidyverse)
require(ggtext)
## read data --------------------------------------------------------------
wdir <- "./chapters/chapter4/"
df <- readxl::read_xlsx(paste0(wdir, "suppl/suppl_2.xlsx"),
                        sheet = "Table S4",
                        skip = 4)

## tidy-up data -----------------------------------------------------------

df$Quant_Idx <-
  factor(df$Quant_Idx,
         levels =  c("GenBank", "RefSeq", "LsRTDv1"))

df$Direction <-
  factor(df$Direction,
         levels =  c("up-regulated", "down-regulated"))

df <-
  df %>%
  group_by(Quant_Idx, Contrasts) %>%
  mutate(label_y = cumsum(N_DE_Gene))


## visualize data (DE genes) -----------------------------------------------

plt1 <- df %>%
  ggplot(aes(x  = Quant_Idx, y  = N_DE_Gene, fill = Direction)) +
  geom_col(position = "stack") +
  geom_text(
    aes(y = label_y, label = N_DE_Gene),
    colour = "black",
    hjust = 1.1,
    angle = 90
  ) +
  labs(x = 'Hours post inoculation',
       y = 'Number of DE genes' ,
       fill = 'Regulation') +
  # scale_x_continuous(breaks = seq(9, 54, 3)) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 13000, 1000),
    limits = c(0, 13000)
  ) +
  scale_fill_manual(values = c(
    "up-regulated" = "#E69F00",
    "down-regulated" = "#56B4E9"
  )) +
  facet_grid( ~ HPI, switch = "both") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = -90),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.5, 0.95),
    legend.direction = "horizontal",
    legend.title = element_blank()
  )

## visualize data (DE transcripts) -----------------------------------------

df <-
  df %>%
  group_by(Quant_Idx, Contrasts) %>%
  mutate(label_y = cumsum(N_DE_Transcript))

plt2 <- df %>%
  ggplot(aes(x  = Quant_Idx, y  = N_DE_Transcript, fill = Direction)) +
  geom_col(position = "stack") +
  geom_text(
    aes(y = label_y, label = N_DE_Transcript),
    colour = "black",
    hjust = 1.1,
    angle = 90
  ) +
  labs(x = 'Hours post inoculation',
       y = 'Number of DE transcripts' ,
       fill = 'Regulation') +
  # scale_x_continuous(breaks = seq(9, 54, 3)) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 13000, 1000),
    limits = c(0, 13000)
  ) +
  scale_fill_manual(values = c(
    "up-regulated" = "#E69F00",
    "down-regulated" = "#56B4E9"
  )) +
  facet_grid( ~ HPI, switch = "both") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = -90),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

# then combine with the top row for final plot
plt <- cowplot::plot_grid(plt1,
                          plt2,
                          labels = "AUTO",
                          label_size = 18,
                          ncol = 1)
# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_DE_genes_tx.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300
)
plt
dev.off()

svg(
  file = paste0(wdir, "figures/time-series_DE_genes_tx.svg"),
  width = 10,
  height = 10,
)
plt
dev.off()
