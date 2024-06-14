rm(list = ls())

library("tidyverse")
library("ComplexHeatmap")

# read dataframe
dat <- readxl::read_xlsx("./suppl/ch4/suppl_1.xlsx",
                         sheet = "Table S96",
                         skip = 4)

df <- dat
CPC2 <- df %>% filter(CPC2 == "*") %>% dplyr::pull(Transcript_ID)
CPAT <- df %>% filter(CPAT == "*") %>% dplyr::pull(Transcript_ID)
FEELnc <-
  df %>% filter(FEELnc == "*") %>% dplyr::pull(Transcript_ID)
No_SwissProt_Hit <-
  df %>% filter(No_SwissProt_Hit == "*") %>% dplyr::pull(Transcript_ID)
TranSuite <- df %>% dplyr::pull(Transcript_ID)

lt = list(
  CPC2 = CPC2,
  CPAT = CPAT,
  FEELnc = FEELnc,
  `No SwissProt Hit` = No_SwissProt_Hit,
  TranSuite = TranSuite
)

list_to_matrix(lt)

m = make_comb_mat(lt)
ss = set_size(m)
cs = comb_size(m)
num.conf.high = sum(comb_size(m, degree = c(5, 4)))
num.conf.mod = sum(comb_size(m, degree = 3))
num.conf.low = sum(comb_size(m, degree = c(2, 1)))

library("circlize")
lgd = Legend(
  at = c(
    paste0("High (", num.conf.high, ")"),
    paste0("Moderate (", num.conf.mod, ")"),
    paste0("Low (", num.conf.low, ")")
  ),
  title = "Conf.\nLevel",
  legend_gp = gpar(fill = c("black", "#5A5A5A", "white")),
  border = "black"
)

ht = UpSet(
  m,
  set_order = order(ss),
  comb_order = order(comb_degree(m), cs),
  pt_size = unit(5, "mm"),
  lwd = 3,
  top_annotation = upset_top_annotation(
    m,
    height = unit(8, "cm"),
    add_numbers = T,
    annotation_name_rot = 90,
    gp = gpar(fill = c(
      "white", "white", "#5A5A5A", "black", "black"
    )[comb_degree(m)])
  ),
  right_annotation = upset_right_annotation(
    m,
    add_numbers = T,
    ylim = c(0, 18000),
    gp = gpar(fill = "black"),
    annotation_name_side = "top",
    axis_param = list(side = "top")
  )
)

svg(
  "figures/ch4/LsRTDv1_lncRNA_tx_upset.svg",
  width = 6,
  height = 4.5
)

draw(ht) %v% draw(lgd, x = unit(0.3, "npc"), y = unit(0.9, "npc"))

dev.off()

summ <-
  df %>%
  group_by(Class_Code) %>%
  summarise(Number = n())

summ <- summ %>%
  mutate(Percentage = round((summ$Number / sum(summ$Number) * 100), digits = 1)) %>%
  mutate(Category = ifelse(
    Class_Code == "i",
    "Intronic",
    ifelse(Class_Code == "x",
           "Antisense",
           "Intergenic")
  )) %>%
  mutate(
    csum = rev(cumsum(rev(Number))),
    pos = Number / 2 + lead(csum, 1),
    pos = if_else(is.na(pos), Number / 2, pos)
  ) %>%
  mutate(
    label = paste(Category, Number, sep = "\n"),
    label = paste0(label, " (", round(Percentage, digits = 1), "%)")
  )

library(ggrepel)
library(ggpattern)

plt <-
  ggplot(summ, aes(
    x = "" ,
    y = Number,
    fill = fct_inorder(Category)
  )) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  geom_text_repel(
    data = summ,
    aes(y = pos, label = paste0(label)),
    nudge_x = 1,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("white", "#5A5A5A", "#D3D3D3")) +
  theme_void(base_size = 16) +
  theme(
    legend.position = "none",
    plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
    plot.title = element_blank()
  )

plt

svg(
  "figures/ch4/LsRTDv1_lncRNA_tx_classes.svg",
  width = 3,
  height = 3
)

plt

dev.off()


df <- dat %>% distinct(Gene_ID, .keep_all = T)

CPC2 <- df %>% filter(CPC2 == "*") %>% dplyr::pull(Transcript_ID)
CPAT <- df %>% filter(CPAT == "*") %>% dplyr::pull(Transcript_ID)
FEELnc <-
  df %>% filter(FEELnc == "*") %>% dplyr::pull(Transcript_ID)
No_SwissProt_Hit <-
  df %>% filter(No_SwissProt_Hit == "*") %>% dplyr::pull(Transcript_ID)
TranSuite <- df %>% dplyr::pull(Transcript_ID)

lt = list(
  CPC2 = CPC2,
  CPAT = CPAT,
  FEELnc = FEELnc,
  `No SwissProt Hit` = No_SwissProt_Hit,
  TranSuite = TranSuite
)

list_to_matrix(lt)

m = make_comb_mat(lt)
ss = set_size(m)
cs = comb_size(m)
num.conf.high = sum(comb_size(m, degree = c(5, 4)))
num.conf.mod = sum(comb_size(m, degree = 3))
num.conf.low = sum(comb_size(m, degree = c(2, 1)))

library("circlize")
lgd = Legend(
  at = c(
    paste0("High (", num.conf.high, ")"),
    paste0("Moderate (", num.conf.mod, ")"),
    paste0("Low (", num.conf.low, ")")
  ),
  title = "Conf.\nLevel",
  legend_gp = gpar(fill = c("black", "#5A5A5A", "white")),
  border = "black"
)

ht = UpSet(
  m,
  set_order = order(ss),
  comb_order = order(comb_degree(m), cs),
  pt_size = unit(5, "mm"),
  lwd = 3,
  top_annotation = upset_top_annotation(
    m,
    height = unit(8, "cm"),
    add_numbers = T,
    annotation_name_rot = 90,
    gp = gpar(fill = c(
      "white", "white", "#5A5A5A", "black", "black"
    )[comb_degree(m)])
  ),
  right_annotation = upset_right_annotation(
    m,
    add_numbers = T,
    ylim = c(0, 18000),
    gp = gpar(fill = "black"),
    annotation_name_side = "top",
    axis_param = list(side = "top")
  )
)

svg(
  "figures/ch4/LsRTDv1_lncRNA_gene_upset.svg",
  width = 8,
  height = 6
)

draw(ht) %v% draw(lgd, x = unit(0.3, "npc"), y = unit(0.9, "npc"))

dev.off()


summ <-
  df %>% distinct(Gene_ID, .keep_all = T) %>%
  group_by(Class_Code) %>%
  summarise(Number = n())

summ <- summ %>%
  mutate(Percentage = round((summ$Number / sum(summ$Number) * 100), digits = 1)) %>%
  mutate(Category = ifelse(
    Class_Code == "i",
    "Intronic",
    ifelse(Class_Code == "x",
           "Antisense",
           "Intergenic")
  )) %>%
  mutate(
    csum = rev(cumsum(rev(Number))),
    pos = Number / 2 + lead(csum, 1),
    pos = if_else(is.na(pos), Number / 2, pos)
  ) %>%
  mutate(
    label = paste(Category, Number, sep = "\n"),
    label = paste0(label, " (", round(Percentage, digits = 1), "%)")
  )

library(ggrepel)

plt <-
  ggplot(summ, aes(
    x = "" ,
    y = Number,
    fill = fct_inorder(Category)
  )) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  geom_text_repel(
    data = summ,
    aes(y = pos, label = paste0(label)),
    size = 5,
    nudge_x = 1,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("white", "#5A5A5A", "#D3D3D3")) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
    plot.title = element_blank()
  )

svg(
  "figures/ch4/LsRTDv1_lncRNA_genes_classes.svg",
  width = 3,
  height = 3
)

plt

dev.off()
