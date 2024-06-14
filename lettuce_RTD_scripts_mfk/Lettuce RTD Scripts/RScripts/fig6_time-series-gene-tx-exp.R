rm(list = ls())
# Gene and Transcript Expression (Lsat_Bcin_timeseries)
# load libraries ----------------------------------------------------------
require(tidyverse)
require(readxl)

# read data --------------------------------------------------------------
wdir <- "./chapters/chapter4/"

df <- readxl::read_xlsx(path = paste0(wdir, "suppl/suppl_2.xlsx"),
                        sheet = "Table S3",
                        skip = 4)
## tidy-up data -----------------------------------------------------------
df <- df %>% select(
  "Quant_Idx",
  "Raw transcripts",
  "Raw genes",
  "Expressed transcripts",
  "Expressed genes"
) %>%
  pivot_longer(!Quant_Idx, names_to = "condition", values_to = "count") %>%
  separate(condition, into = c("exp", "level"), sep = " (?=[^ ]+$)") %>%
  mutate(exp = ifelse(exp == "Raw", "Total", "Expressed"))

df$Quant_Idx <-
  factor(df$Quant_Idx,
         levels =  c("GenBank", "RefSeq", "LsRTDv1"))

labs <- c("Gene Expression", "Transcript Expression")
names(labs) <- c("genes", "transcripts")

## visualize data --------------------------------------------------------------
plt <- df %>%
  ggplot(aes(x = Quant_Idx, y = count, fill = exp)) +
  geom_col(color = "black", position = 'identity') +
  geom_text(aes(label = count, colour = exp),
            vjust = 1.6) +
  scale_fill_manual(values = c("black", "white")) +
  scale_colour_manual(values = c("white", "black")) +
  scale_y_continuous(
    limits = c(0, 200000),
    breaks = seq(0, 200000, by = 50000),
    expand = c(0, 0)
  ) +
  facet_grid(. ~ level, labeller = labeller(level = labs)) +
  labs(x = 'Transcriptome', y = 'Count number') +
  theme_classic(base_size = 16) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    strip.text.x = element_text(color = "black",
                                face = "bold")
  )

# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_gene_tx_exp.png"),
  width = 6,
  height = 4.5,
  units = "in",
  res = 300
)

plt

dev.off()

svg(
  file = paste0(wdir, "figures/time-series_gene_tx_exp.svg"),
  width = 6,
  height = 4.5
)

plt

dev.off()