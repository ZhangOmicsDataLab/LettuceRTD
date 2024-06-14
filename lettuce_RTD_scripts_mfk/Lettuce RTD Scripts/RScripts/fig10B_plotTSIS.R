rm(list = ls())
# Plot expression profile of time-series isoform switches
# plot function -----------------------------------------------------------
plotTSIS <-
  function(Gene.TPM,
           Transcript.TPM,
           scores = NULL,
           func.annot = NULL,
           GeneID = NULL,
           prot.feature = NULL,
           t.points.cutoff = NULL,
           pval.cutoff = NULL,
           cor.cutoff = NULL,
           prob.cutoff = 0.5,
           Tx = NULL) {
    # load libraries ----------------------------------------------------------
    require(tidyverse)
    require(xlsx)
    require(ggtext)
    # read R function ---------------------------------------------------------
    wdir <- "./chapters/chapter4/"
    source(paste0(wdir, "scripts/func/read_DDD_results.R"))
    
    # tidy-up TSIS scores -----------------------------------------------------
    scores <-
      scores |>
      filter(treatment == Tx) |>
      filter(
        before.t.points >= t.points.cutoff &
          after.t.points >= t.points.cutoff &
          before.pval < pval.cutoff &
          after.pval < pval.cutoff &
          abs(cor) >= cor.cutoff
      ) |>
      mutate(gene = gsub("\\.[0-9]*$", "", iso1)) |>
      relocate(gene, .before = iso1) |>
      filter(gene == GeneID)
    
    # pull isoform ids --------------------------------------------------------
    iso1 = scores |> pull(iso1)
    iso2 = scores |> pull(iso2)
    iso = c(iso1, iso2) %>% unique()
    
    # create title and subtitle for the plot ---------------------------------
    if (GeneID %in% func.annot$gene.id) {
      func.annot <- func.annot |>
        distinct(gene.id, .keep_all = T) |>
        filter(gene.id == GeneID)
      
      if (is.na(func.annot$ahrd.desc)) {
        title = paste0(func.annot$gene.id)
        subtitle = NULL
      } else {
        title = paste0(func.annot$gene.id)
        subtitle = paste0(
          func.annot$ahrd.desc
          # "\n",
          # "BLAST Hit: ",
          # func.annot$source,
          # "|",
          # func.annot$organism,
          # "|",
          # func.annot$primary.acc,
          # " [",
          # func.annot$ahrd.quality.code,
          # "]"
        )
      }
    } else{
      title = paste0(GeneID)
      subtitle = NULL
    }
    
    print(paste("Gene ID:", title))
    print(paste("Gene description:", subtitle))
    print(paste(c("Isoform IDs:", iso), collapse = " "))
    
    # transcript feature ------------------------------------------------------
    prot.feature <- prot.feature |>
      filter(Transcript_ID %in% iso) |>
      select(Transcript_ID, Coding_potentiality, Features, Translation) |>
      group_by(Transcript_ID) |>
      mutate(legend.label = paste0(
        Transcript_ID,
        " (",
        Coding_potentiality,
        if (Features != "-") {
          paste0(" [", Features, "]")
        },
        ", ",
        paste0(nchar(Translation), "aa"),
        ")"
      )) |>
      rename(id = Transcript_ID)
    # read expression files ---------------------------------------------------
    Gene.TPM <- Gene.TPM |>
      rename(geneid = X) |>
      filter(geneid == GeneID) |>
      pivot_longer(!geneid, names_to = "sample", values_to = "TPM")  |>
      separate(sample, c("treatment", "hpi", "brep"), sep = "\\.")  |>
      filter(treatment == paste0(Tx)) |>
      rename(id = "geneid")
    
    Transcript.TPM <-
      Transcript.TPM  |>
      rename(txid = X) |>
      filter(txid %in% iso) |>
      pivot_longer(!txid, names_to = "sample", values_to = "TPM")  |>
      separate(sample, c("treatment", "hpi", "brep"), sep = "\\.")  |>
      filter(treatment == paste0(Tx)) |>
      rename(id = "txid")
    
    Gene.Tx.TPM <-
      rbind(Gene.TPM, Transcript.TPM)  |>
      mutate(hpi = as.numeric(gsub('[[:alpha:]]', '', hpi)))  |>
      arrange(desc(hpi)) |>
      group_by(id, treatment, hpi)  |>
      mutate(
        mean.TPM = mean(TPM),
        SE.TPM = sd(TPM) /
          sqrt(3),
        high.TPM = mean.TPM + SE.TPM,
        low.TPM = mean.TPM - SE.TPM
      )
    
    
    # color schema ------------------------------------------------------------
    colours <-
      c(
        "#009E73",
        "#CC79A7",
        "#D55E00",
        "#0072B2",
        "#F0E442",
        "#E69F00",
        "#56B4E9",
        "#D3D3D3",
        "#5A5A5A"
      )
    # plot dataframe ----------------------------------------------------------
    
    Gene.Tx.TPM %>% ggplot2::ggplot() +
      geom_line(aes(x = hpi, y = mean.TPM, colour = id), linewidth = 1) +
      geom_point(aes(x = hpi, y = TPM, colour = id), size = 1) +
      geom_ribbon(aes(
        x = hpi,
        ymin = low.TPM,
        ymax = high.TPM,
        fill = id
      ), alpha = 0.3) +
      geom_point(data = scores, aes(x = x.value, y = y.value), size = 4, shape = 10, stroke = 2) +
      scale_fill_manual(
        values = c("#5A5A5A", colours[1:length(iso)], 'black'),
        breaks = c(paste0(GeneID), paste0(iso), 'switch_points'),
        labels = c(
          'Total gene',
          gsub(".*\\.", ".", prot.feature$legend.label),
          'switch_points'
        )
      ) +
      scale_color_manual(
        values = c("#5A5A5A", colours[1:length(iso)], 'black'),
        breaks = c(paste0(GeneID), paste0(iso), 'switch_points'),
        labels = c(
          'Total gene',
          gsub(".*\\.", ".", prot.feature$legend.label),
          'switch_points'
        )
      ) +
      # annotate(
      #   'text',
      #   x = 4 + scores$x.value,
      #   y = scores$y.value + max(Gene.Tx.TPM$mean.TPM) / 20,
      #   label = paste0(
      #     'prob=',
      #     round(scores$prob, 2),
      #     '\n',
      #     'diff=',
      #     round(scores$diff, 1),
      #     '\n',
      #     'cor=',
      #     round(scores$cor, 2)
      #   )
      # ) +
      scale_x_continuous(
        name = "Hours post inoculation",
        breaks = unique(Gene.Tx.TPM$hpi),
        labels = unique(Gene.Tx.TPM$hpi)) +
      scale_y_continuous(
        name = "TPM"
      ) +
      labs(title = title, subtitle = subtitle) +
      theme_classic(base_size = 18) +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.text = element_markdown()
      )
  }

wdir <- "./chapters/chapter4/"
# read data required ------------------------------------------------------
df1 <- readxl::read_xlsx(path = paste0(wdir, "suppl/Table S20.xlsx"),
                         skip = 4)

df2 <-
  read.csv(
    "U:/Lettuce_RTD_analysis/P08E96_references/LettuceRTDv1/func_annot/AHRD/AHRD_func_annot.csv"
  )

# set parameters ----------------------------------------------------------
root.dir = "U:/Lettuce_RTD_analysis" # root files directory
df <- read.csv(paste(root.dir, "3D_output_dirs.csv", sep = "/"))
dataset <- unique(df$dataset)
ref.idx <- unique(df$ref.idx)
# read R function ---------------------------------------------------------
source("./analyses/ch4/func/read_DDD_results.R")
# read TPM expression files -----------------------------------------------

df3 <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "Gene TPM.csv"
)

df4 <- read.DDD.results(
  zip.path = df$zip.path[3],
  job.id = df$job.ids[3],
  file.name = "Transcript TPM.csv"
)

# read feature table ------------------------------------------------------

df5 <- readxl::read_xlsx(path = paste0(wdir, "suppl/Table S8.xlsx"),
                         skip = 4)

gene.list <- c(
  "LSATv11_C01_002961", # * may involve in cell-wall organisation and defence response
  "LSATv11_C01_003796", # *
  "LSATv11_C01_004406",
  "LSATv11_C01_008212", # * transport
  "LSATv11_C01_008490", # * fatty acid alpha-oxidation, pyruvate decarboxylase, lyase activity, peroxisome
  "LSATv11_C02_011992",
  "LSATv11_C02_012905",
  "LSATv11_C02_013302", # oxidoreductase activity
  "LSATv11_C02_013512", # transferase activity
  "LSATv11_C02_014177", # * RNA methyltransferase activity
  "LSATv11_C02_014816", # * zinc ion binding
  "LSATv11_C03_015900", # oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor
  "LSATv11_C03_016574",
  "LSATv11_C03_018385", # amino acid transport, gamma-aminobutyric acid transport
  "LSATv11_C03_018427", # protein phosphorylation
  "LSATv11_C03_018857", # *** autophagosome assembly, autophagy of mitochondrion, plant defense, stress response
  "LSATv11_C03_018992", # xenobiotic transmembrane transporter activity
  "LSATv11_C03_019896", # ** protein phosphorylation, MAPKKK YODA-like, 	MAPK cascade, defense response to fungus
  "LSATv11_C03_020327", # * Putative actin cross-linking
  "LSATv11_C04_028665", # * MAPK cascade, protein phosphorylation, MAPKKK YODA-like
  "LSATv11_C04_031355", # protein phosphorylation, MAPK cascade
  "LSATv11_C05_035728", # protein dephosphorylation
  "LSATv11_C05_037687", # translation, structural constituent of ribosome
  "LSATv11_C05_038031", # ***
  "LSATv11_C05_040232", # * protein phosphorylation
  "LSATv11_C06_041084",
  "LSATv11_C06_041552",
  "LSATv11_C07_046210",
  "LSATv11_C07_048500", # positive regulation of ATP-dependent activity, protein folding
  "LSATv11_C07_049829", # pyridoxal 5'-phosphate salvage
  "LSATv11_C07_051157", 
  "LSATv11_C08_052712",
  "LSATv11_C08_052726",
  "LSATv11_C08_053557", # ** clathrin coat assembly, clathrin-dependent endocytosis
  "LSATv11_C08_054192",
  "LSATv11_C08_055033",
  "LSATv11_C08_059159", # * histone H3-K27 trimethylation, heterochromatin formation
  "LSATv11_C09_059900", # *** auxin-activated signaling pathway, 	regulation of DNA-templated transcription
  "LSATv11_C09_063822", 
  "LSATv11_C09_063875"
)

for (i in 1:length(gene.list)) {
  png(
    file = paste0("./figures/ch4/iso-switch_plots/",
                  gene.list[i],
                  ".png"),
    units = "in",
    width = 8,
    height = 6,
    res = 300
  )
  print(
    plotTSIS(
      Gene.TPM = df3,
      Transcript.TPM = df4,
      scores = df1,
      func.annot = df2,
      GeneID = gene.list[i],
      prot.feature = df5,
      t.points.cutoff = 2,
      pval.cutoff = 0.01,
      cor.cutoff = 0,
      prob.cutoff = 0.5,
      Tx = "INOC"
    )
  )
  dev.off()
}

gene.list <- c(
  "LSATv11_C01_002961", # *** may involve in cell-wall organisation and defence response, carbohydrate metabolic process
  "LSATv11_C03_018857", # *** autophagosome assembly, autophagy of mitochondrion, plant defense, stress response, autophagy
  "LSATv11_C05_038031", # ***
  "LSATv11_C08_053557", # ** clathrin coat assembly, clathrin-dependent endocytosis, photosynthetic electron transport in photosystem I
  "LSATv11_C08_059159", # * histone H3-K27 trimethylation, heterochromatin formation, histone modification, chromatin remodeling
  "LSATv11_C09_059900"  # *** auxin-activated signaling pathway, 	regulation of DNA-templated transcription
)

p1 <- plotTSIS(
  Gene.TPM = df3,
  Transcript.TPM = df4,
  scores = df1,
  func.annot = df2,
  GeneID = gene.list[1],
  prot.feature = df5,
  t.points.cutoff = 2,
  pval.cutoff = 0.01,
  cor.cutoff = 0,
  prob.cutoff = 0.5,
  Tx = "INOC"
)

p2 <- plotTSIS(
  Gene.TPM = df3,
  Transcript.TPM = df4,
  scores = df1,
  func.annot = df2,
  GeneID = gene.list[2],
  prot.feature = df5,
  t.points.cutoff = 2,
  pval.cutoff = 0.01,
  cor.cutoff = 0,
  prob.cutoff = 0.5,
  Tx = "INOC"
)
p3 <- plotTSIS(
  Gene.TPM = df3,
  Transcript.TPM = df4,
  scores = df1,
  func.annot = df2,
  GeneID = gene.list[3],
  prot.feature = df5,
  t.points.cutoff = 2,
  pval.cutoff = 0.01,
  cor.cutoff = 0,
  prob.cutoff = 0.5,
  Tx = "INOC"
)
p4 <- plotTSIS(
  Gene.TPM = df3,
  Transcript.TPM = df4,
  scores = df1,
  func.annot = df2,
  GeneID = gene.list[4],
  prot.feature = df5,
  t.points.cutoff = 2,
  pval.cutoff = 0.01,
  cor.cutoff = 0,
  prob.cutoff = 0.5,
  Tx = "INOC"
)
p5 <- plotTSIS(
  Gene.TPM = df3,
  Transcript.TPM = df4,
  scores = df1,
  func.annot = df2,
  GeneID = gene.list[5],
  prot.feature = df5,
  t.points.cutoff = 2,
  pval.cutoff = 0.01,
  cor.cutoff = 0,
  prob.cutoff = 0.5,
  Tx = "INOC"
)
p6 <- plotTSIS(
  Gene.TPM = df3,
  Transcript.TPM = df4,
  scores = df1,
  func.annot = df2,
  GeneID = gene.list[6],
  prot.feature = df5,
  t.points.cutoff = 2,
  pval.cutoff = 0.01,
  cor.cutoff = 0,
  prob.cutoff = 0.5,
  Tx = "INOC"
)

# save plot ---------------------------------------------------------------
png(
  file = paste0(wdir, "figures/iso_switch_exp_vert.png"),
  width = 12,
  height = 16,
  units = "in",
  res = 300
)

cowplot::plot_grid(
  p1,
  p2,
  p3,
  p4,
  p5,
  p6,
  nrow = 3,
  ncol = 2)

dev.off()

svg(
  file = paste0(wdir, "figures/iso_switch_exp_vert.svg"),
  width = 12,
  height = 16
  )

cowplot::plot_grid(
  p1,
  p2,
  p3,
  p4,
  p5,
  p6,
  nrow = 3,
  ncol = 2)

dev.off()

png(
  file = paste0(wdir, "figures/iso_switch_exp_horiz.png"),
  width = 16,
  height = 10,
  units = "in",
  res = 300
)

cowplot::plot_grid(p1,
                   p2,
                   p5,
                   p6,
                   nrow = 2,
                   ncol = 2)

dev.off()

svg(
  file = paste0(wdir, "figures/iso_switch_exp_horiz.svg"),
  width = 16,
  height = 10
)

cowplot::plot_grid(p1,
                   p2,
                   p5,
                   p6,
                   nrow = 2,
                   ncol = 2)

dev.off()
