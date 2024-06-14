rm(list = ls())
# Mapping Efficiency (Lsat_DivSet_Sscl)
require(tidyverse)
require(ggbeeswarm)
require(rstatix)
# read data --------------------------------------------------------------
wdir <- "./Chapters/ch4/"
metadata <- readxl::read_xlsx(
  path = paste0(wdir, "suppl/suppl_3.xlsx"),
  sheet = "Table S1",
  skip = 4
)
multiqc_salmon <- readxl::read_xlsx(
  path = paste0(wdir, "suppl/suppl_3.xlsx"),
  sheet = "Table S2",
  skip = 4
)

## tidy-up input and rename samples ----------------------------------------
multiqc_salmon <- multiqc_salmon |>
  select(sample, percent_mapped, ref_tx_idx) |>
  mutate(
    ref_tx_idx = ifelse(
      ref_tx_idx == "Quant_GCA_002870075.4_Lsat_Salinas_v11_GCA_001857865.1_ASM185786v1_merged_salmon_index",
      "GenBank",
      ifelse(
        ref_tx_idx == "Quant_GCF_002870075.4_Lsat_Salinas_v11_GCA_001857865.1_ASM185786v1_merged_salmon_index",
        "RefSeq",
        "LsRTDv1"
      )
    )
  )

df <-
  readxl::read_xlsx(path = paste0(wdir, "suppl/metadata_simplified.xlsx"),
                    range = "A1:D74")
sample.excl = df |>
  filter(`is.SRA.avail?` == FALSE) |>
  pull(sample)

multiqc_salmon <- multiqc_salmon |>
  filter(!sample %in% sample.excl)


str(multiqc_salmon)
summary(multiqc_salmon)
# Check the number of paired observations
xtabs( ~ sample + ref_tx_idx, multiqc_salmon)
# Histogram of difference data
GenBank = multiqc_salmon$percent_mapped[multiqc_salmon$ref_tx_idx == "GenBank"]
RefSeq = multiqc_salmon$percent_mapped[multiqc_salmon$ref_tx_idx == "RefSeq"]
LsRTDv1 = multiqc_salmon$percent_mapped[multiqc_salmon$ref_tx_idx == "LsRTDv1"]
LsRTDv1vsGenBank.Diff = LsRTDv1 - GenBank
LsRTDv1vsRefSeq.Diff = LsRTDv1 - RefSeq

library(rcompanion)
library(lsr)
plotNormalHistogram(LsRTDv1vsGenBank.Diff,
                    xlab = "Difference (LsRTDv1 - GenBank)")
plotNormalHistogram(LsRTDv1vsRefSeq.Diff,
                    xlab = "Difference (LsRTDv1 - RefSeq)")
# Plot the paired data
GenBank = multiqc_salmon$percent_mapped[multiqc_salmon$ref_tx_idx == "GenBank"]
RefSeq = multiqc_salmon$percent_mapped[multiqc_salmon$ref_tx_idx == "RefSeq"]
LsRTDv1 = multiqc_salmon$percent_mapped[multiqc_salmon$ref_tx_idx == "LsRTDv1"]
Samples = multiqc_salmon$sample[multiqc_salmon$ref_tx_idx == "LsRTDv1"]

plot(
  GenBank,
  jitter(LsRTDv1),
  # jitter offsets points so you can see them all
  pch = 16,
  # shape of points
  cex = 1.0,
  # size of points
  xlim = c(50, 110),
  # limits of x-axis
  ylim = c(50, 110),
  # limits of y-axis
  xlab = "GenBank",
  # label for x-axis
  ylab = "LsRTDv1"
)             # label for y-axis

text(GenBank,
     LsRTDv1,
     labels = Samples,
     # Label location and text
     
     pos = 3,
     cex = 1.0)               # Label text position and size

abline(0, 1, col = "blue", lwd = 2)     # line with intercept of 0 and slope of 1

barplot(
  LsRTDv1vsGenBank.Diff,
  # variable to plot
  col = "dark gray",
  # color of bars
  xlab = "Samples",
  # x-axis label
  ylab = "Difference (LsRTDv1 - GenBank)",
  # y-axis label
  names.arg = Samples
)                        # labels for bars

plot(
  RefSeq,
  jitter(LsRTDv1),
  # jitter offsets points so you can see them all
  pch = 16,
  # shape of points
  cex = 1.0,
  # size of points
  xlim = c(50, 110),
  # limits of x-axis
  ylim = c(50, 110),
  # limits of y-axis
  xlab = "RefSeq",
  # label for x-axis
  ylab = "LsRTDv1"
)           # label for y-axis

text(RefSeq,
     LsRTDv1,
     labels = Samples,
     # Label location and text
     
     pos = 3,
     cex = 1.0)               # Label text position and size

abline(0, 1, col = "blue", lwd = 2)     # line with intercept of 0 and slope of 1

barplot(
  LsRTDv1vsRefSeq.Diff,
  # variable to plot
  col = "dark gray",
  # color of bars
  xlab = "Samples",
  # x-axis label
  ylab = "Difference (LsRTDv1 - GenBank)",
  # y-axis label
  names.arg = Samples
)                        # labels for bars


df <- multiqc_salmon |>
  subset(ref_tx_idx != "GenBank") |>
  mutate(ref_tx_idx = as.factor(as.character(ref_tx_idx)))

str(df)
t.test(percent_mapped ~ ref_tx_idx,
       data   = df,
       paired = TRUE)

cohensD(percent_mapped ~ ref_tx_idx,
        data   = df,
        method = "paired")

df <- multiqc_salmon |>
  subset(ref_tx_idx != "RefSeq") |>
  mutate(ref_tx_idx = as.factor(as.character(ref_tx_idx)))

str(df)
t.test(percent_mapped ~ ref_tx_idx,
       data   = df,
       paired = TRUE)

cohensD(percent_mapped ~ ref_tx_idx,
        data   = df,
        method = "paired")
# Convert id and time into factor variables
multiqc_salmon <- multiqc_salmon %>%
  convert_as_factor(sample, ref_tx_idx)
# Inspect some random rows of the data by groups
set.seed(123)
multiqc_salmon %>% sample_n_by(ref_tx_idx, size = 2)

# Pairwise comparisons between time points at each group levels
# Paired t-test is used because we have repeated measures by time
library(rstatix)
library(ggpubr)
stat.test <- multiqc_salmon %>%
  pairwise_t_test(percent_mapped ~ ref_tx_idx,
                  paired = TRUE,
                  p.adjust.method = "bonferroni") %>%
  select(-df,-statistic,-p) # Remove details
stat.test

# Create the plot
bxp <- ggboxplot(
  multiqc_salmon,
  x = "ref_tx_idx",
  y = "percent_mapped",
  color = "ref_tx_idx",
  palette = "jco"
)
# Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = "ref_tx_idx")
bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif",
                         step.increase = 0.08)
# Hide ns
bxp + stat_pvalue_manual(
  stat.test,
  label = "p.adj.signif",
  step.increase = 0.08,
  hide.ns = TRUE,
  tip.length = 0
)


### relevel values to make plot ordered
multiqc_salmon$ref_tx_idx <-
  factor(multiqc_salmon$ref_tx_idx,
         # Change ordering manually
         levels = c("GenBank", "RefSeq", "LsRTDv1"))

## plot mapping efficiency -------------------------------------------------
library(superb)
plt <- multiqc_salmon |>
  ggplot2::ggplot(aes(x = ref_tx_idx, y = percent_mapped)) +
  # geom_boxplot(width = 0.2, outlier.color = NA) +
  gghalves::geom_half_point(
    range_scale = .75,
    size = 2,
    alpha = 0.3,
    transformation = position_quasirandom(method = "quasirandom"),
    position = position_nudge(x = .01),
  ) +
  gghalves::geom_half_boxplot(outlier.color = NA,
                              position = position_nudge(x = -.01)) +
  showSignificance(c(2, 3), 78, 1, "**** \n Cohen's d = 1.89",
                   textParams    = list(size = 4)) +
  showSignificance(c(1, 3), 72, 1, "**** \n Cohen's d = 5.22", textParams    = list(size = 4)) +
  scale_y_continuous(
    limits = c(print(
      min(multiqc_salmon$percent_mapped) %/% 10 * 10
    ), 100),
    breaks = seq(print(
      min(multiqc_salmon$percent_mapped) %/% 10 * 10
    ), 100, by = 10),
    expand = c(0, 0)
  ) +
  labs(y = "Mapped reads (%)", x = "Transcriptome") +
  theme_classic(base_size = 16)
plt

# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/divset_mapping_rate.png"),
  width = 6,
  height = 4.5,
  units = "in",
  res = 300
)
plt
dev.off()

svg(
  file = paste0(wdir, "figures/divset_mapping_rate.svg"),
  width = 6,
  height = 4.5,
)
plt
dev.off()
