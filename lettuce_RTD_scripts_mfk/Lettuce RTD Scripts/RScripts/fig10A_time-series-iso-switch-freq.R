rm(list = ls())
# Frequency of significant isoform switches
# over the course of B. cinerea infection.
# load library ------------------------------------------------------------
library(tidyverse)
library(ggtext)
wdir <- "./chapters/chapter4/"
# read dataframe ----------------------------------------------------------
df <- readxl::read_xlsx(
  path = paste0(wdir, "suppl/suppl_2.xlsx"),
                      sheet = "Table S8",
                      skip = 4)
# tidy-up dataframe -------------------------------------------------------
t.points.cutoff = 2
pval.cutoff = 0.01
cor.cutoff = 0

df <- df |>
  mutate(gene.id = gsub("\\..*", "", iso1)) |> 
  filter(
    before.t.points >= t.points.cutoff &
      after.t.points >= t.points.cutoff &
      before.pval < pval.cutoff &
      after.pval < pval.cutoff &
      abs(cor) > cor.cutoff
  )

# visualize --------------------------------------------------------------
plt <- ggplot(data = df, aes(x = x.value))  +
  geom_histogram(
    aes(fill = treatment),
    bins = 14,
    position = "dodge",
    alpha = 0.7,
    colour = "black"
  ) +
  stat_bin(
    bins = 14, geom = "text", color='black',
    aes(label=after_stat(ifelse(count == 0, "", count)), group=treatment), position=position_dodge(),vjust = -0.5)+
  scale_x_continuous(
    name = "Hours post inoculation",
    limits = c(9, 48),
    breaks = seq(9, 48, by = 3),
    expand = c(0, 1),
    oob = scales::squish_infinite
  ) +
  scale_y_continuous(name = "Number of isoform switches",
                     limits = c(0, 70),
                     breaks = seq(0, 70, by = 10),
                     expand = c(0, 0)) +
  scale_fill_manual(
    values = c("MOCK" = "#009E73", "INOC" = "#CC79A7"),
    labels = c("MOCK" = "Mock", "INOC" = "*B. cinerea*")
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.7, 0.9),
    legend.direction = "horizontal",
    legend.text = element_markdown()
  )

plt
# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/iso_switch_freq.png"),
  width = 6,
  height = 4.5,
  units = "in",
  res = 300
)

plt

dev.off()

svg(
  file = paste0(wdir, "figures/iso_switch_freq.svg"),
  width = 6,
  height = 4.5
)

plt

dev.off()
