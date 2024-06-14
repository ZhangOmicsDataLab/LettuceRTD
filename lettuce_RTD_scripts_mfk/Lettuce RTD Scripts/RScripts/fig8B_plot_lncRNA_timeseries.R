plot.timeseries <- function(gene.id) {
  library(ggplot2)
  library(ggtext)
  ## read R function --------------------------------------------------------
  source("./analyses/ch4/func/read_DDD_results.R")
  ## set directories ---------------------------------------------------------
  root.dir = "U:/Lettuce_RTD_analysis" # root files directory
  output.dirs <-
    read_csv(paste(root.dir, "3D_output_dirs.csv", sep = "/"),
             show_col_types = F) # 3D RNA-seq output directories
  # gene expression matrix --------------------------------------------------
  Gene.TPM <- read.DDD.results(
    zip.path = output.dirs$zip.path[3],
    job.id = output.dirs$job.ids[3],
    file.name = "Gene TPM.csv"
  )
  # tidy-up gene expression matrix ------------------------------------------
  exp.target.gene <- Gene.TPM |>
    filter(X == {{ gene.id }}) |>
    pivot_longer(!X, names_to = "contrast", values_to = "TPM") |>
    separate(contrast,
             into = c("treatment", "HPI", "brep"),
             sep = "\\.")
  # Remove 'T' and convert to numeric
  exp.target.gene$HPI <-
    as.numeric(gsub("T", "", exp.target.gene$HPI))
  # summarize expression data
  sum.exp.target.gene <- exp.target.gene |>
    Rmisc::summarySE(measurevar = "TPM",
                     groupvars = c("X", "treatment", "HPI")) |>
    mutate(y.high = TPM + se,
           y.low = TPM - se)
  
  #  create title and subtitle for the plot ---------------------------------
    title = paste0(gene.id)
    subtitle = NULL
  # visualize gene expression -----------------------------------------------
  exp.target.gene |>
    ggplot(aes(x = HPI, y = TPM, colour = treatment)) +
    geom_line(data = sum.exp.target.gene,
              mapping = aes(x = HPI, y = TPM, colour = treatment)) +
    geom_ribbon(
      data = sum.exp.target.gene,
      mapping = aes(
        x = HPI,
        ymin = y.low,
        ymax = y.high,
        fill = treatment
      ),
      alpha = 0.3,
      colour = NA
    ) +
    geom_point(alpha = 0.7) +
    scale_color_manual(
      values = c("CONTROL" = "#009E73", "INOC" = "#CC79A7"),
      labels = c("CONTROL" = "Mock", "INOC" = "*B. cinerea*")
    ) +
    scale_fill_manual(
      values = c("CONTROL" = "#009E73", "INOC" = "#CC79A7"),
      labels = c("CONTROL" = "Mock", "INOC" = "*B. cinerea*")
    ) +
    scale_x_continuous(
      name = "Hours post inoculation",
      breaks = seq(9, 48, by = 3)) +
    labs(title = title, subtitle = subtitle) +
    theme_classic(base_size = 18) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_markdown()
          )
}
