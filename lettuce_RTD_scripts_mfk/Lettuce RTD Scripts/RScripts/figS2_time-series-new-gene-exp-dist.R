# Distribution of average expression of new vs old genes
require(tidyverse)
require(ggtext)
# read R function ---------------------------------------------------------
wdir <- "./chapters/chapter4/"
source(paste0(wdir, "scripts/func/read_DDD_results.R"))

## set parameters ---------------------------------------------------------
root.dir = "U:/Lettuce_RTD_analysis" # root files directory
df <-
  read_csv(paste(root.dir, "3D_output_dirs.csv", sep = "/")) # 3D RNA-seq output directories

## extract expressed genes ids --------------------------------------------
DE_gene_testing_stats <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "DE gene testing statistics.csv"
)
target.high <- unique(DE_gene_testing_stats$target)
rm(DE_gene_testing_stats)
## gene expression matrix --------------------------------------------------
Gene.TPM <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "Gene TPM.csv"
)

## new gene ids ------------------------------------------------------------
new.genes <- df <- readxl::read_xlsx(
  path = paste0(wdir, "suppl/Table S7.xlsx"),
  skip = 4)
new.genes <- new.genes$Gene_ID

## tidy-up gene expression matrix ------------------------------------------
Gene.TPM.target.high <- Gene.TPM |>
  filter(X %in% target.high) |>
  pivot_longer(!X, names_to = "contrast", values_to = "TPM") |>
  separate(contrast,
           into = c("treatment", "HPI", "brep"),
           sep = "\\.") |>
  group_by(X, treatment) |>
  summarise(mean.TPM = mean(TPM)) |>
  mutate(NEW = ifelse(X %in% new.genes, "TRUE", "FALSE"))
## visualize data ------------------------------------------------------------
plt <- ggplot() +
  gghalves::geom_half_violin(
    data = Gene.TPM.target.high |> filter(NEW == "TRUE"),
    mapping = aes(
      x = treatment,
      y = log2(mean.TPM+0.001),
      fill = NEW,
      color = NEW
    ),
    position = position_nudge(x = -.01),
    side = "l"
  ) +
  gghalves::geom_half_violin(
    data = Gene.TPM.target.high |> filter(NEW == "FALSE"),
    mapping = aes(
      x = treatment,
      y = log2(mean.TPM+0.001),
      fill = NEW,
      color = NEW
    ),
    position = position_nudge(x = .01),
    side = "r"
  ) +
  gghalves::geom_half_boxplot(
    data = Gene.TPM.target.high |> filter(NEW == "TRUE"),
    mapping = aes(
      x = treatment,
      y = log2(mean.TPM+0.001),
      fill = NEW
    ),
    width = .2,
    position = position_nudge(x = -.01),
    side = "l",
    outlier.shape = NA,
    # remove outlier points
    errorbar.draw = F,
    lwd = 1 # change thickness of box line
  ) +
  gghalves::geom_half_boxplot(
    data = Gene.TPM.target.high |> filter(NEW == "FALSE"),
    mapping = aes(
      x = treatment,
      y = log2(mean.TPM+0.001),
      fill = NEW
    ),
    width = .2,
    position = position_nudge(x = .01),
    side = "r",
    outlier.shape = NA,
    # remove outlier points
    errorbar.draw = F,
    lwd = 1 # change thickness of box line
  ) +
  scale_fill_manual(
    values = c("TRUE" = "grey80", "FALSE" = "grey50"),
    labels = c("TRUE" = "Novel", "FALSE" = "Pre-existing")
  ) +
  scale_color_manual(
    values = c("TRUE" = "grey80", "FALSE" = "grey50"),
    labels = c("TRUE" = "Novel", "FALSE" = "Pre-existing")
  ) +
  scale_x_discrete("Treatment", labels = c("CONTROL" = "Mock", "INOC" = "*B. cinerea*")) +
  scale_y_continuous("log2(TPM + 0.001)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_markdown())
plt


# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/time-series_new_gene_exp_dist.png"),
  width = 6,
  height = 4.5,
  units = "in",
  res = 300
)

plt

dev.off()

svg(
  file = paste0(wdir, "figures/time-series_new_gene_exp_dist.svg"),
  width = 6,
  height = 4.5
)

plt

dev.off()
