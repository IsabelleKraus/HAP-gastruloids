library(ggplot2)
library(dplyr)
library(openxlsx)
source("scripts/OMG_colors.r")

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hap   <- read.delim(file.path(sd_dir, "HAP_UMAP_expression.tsv.gz"), check.names = FALSE)
ko    <- read.delim(file.path(sd_dir, "HIF1AKO_on_HAP_UMAP_expression.tsv.gz"), check.names = FALSE)
query <- read.delim(file.path(sd_dir, "OMG_query_UMAP.tsv.gz"))

# numbered labels, e.g. "Midbrain (3)"; factor levels fix the legend order
add_num <- function(df) {
  df$cell_type_num <- factor(unname(celltype_order[df$cell_state]),
                             levels = unname(celltype_order))
  df
}
hap <- add_num(hap)
ko  <- add_num(ko)

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# UMAP coloured and numbered by OMG cell state
plot_states <- function(df, bg = NULL) {
  centroids <- df %>% filter(!is.na(cell_type_num)) %>%
    group_by(cell_type_num, state_id) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
  p <- ggplot()
  if (!is.null(bg)) p <- p + geom_point(data = bg, aes(umap_1, umap_2),
                                        colour = "grey82", size = 3, alpha = 0.5)
  p +
    geom_point(data = df, aes(umap_1, umap_2, colour = cell_type_num), size = 3) +
    scale_colour_manual(name = "", values = cell_type_colored_numbered, na.value = "gray80") +
    geom_text(data = centroids, aes(umap_1, umap_2, label = state_id),
              size = 6, fontface = "bold", colour = "black") +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    coord_fixed() + theme_umap
}


## Fig 5E: HAP WT UMAP 
ggsave(file.path(out_dir, "Fig5E_HAP_WT_UMAP.pdf"), plot_states(hap), width = 10, height = 7)

## Fig 5E: HAP Hif1a KO projected onto the HAP UMAP (06: UMAP_HIF1AKO_on_HAP.pdf)
ggsave(file.path(out_dir, "Fig5E_HIF1AKO_on_HAP_UMAP.pdf"), plot_states(ko, bg = hap),
       width = 10, height = 7)

## Fig 5E: composition bar WT vs Hif1a KO 
##   percentage of cells per predicted OMG cell state within each condition
comp_table <- function(conds) {
  query %>%
    filter(condition %in% conds) %>%
    count(condition, cell_state, name = "count") %>%
    group_by(condition) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
}
plot_composition <- function(conds, labels = conds) {
  comp <- comp_table(conds)
  present <- intersect(names(cell_type_colored), unique(comp$cell_state))   # legend order as in 04
  comp$cell_type_num <- factor(unname(celltype_order[comp$cell_state]),
                               levels = unname(celltype_order[present]))
  comp$condition <- factor(comp$condition, levels = conds, labels = labels)
  ggplot(comp, aes(condition, percentage, fill = cell_type_num)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cell_type_colored_numbered) +
    labs(x = "", y = "Percentage of cells", fill = "Cell states") +
    theme_classic() +
    theme(axis.text.x  = element_text(size = 16, angle = 45, hjust = 1, colour = "black"),
          axis.text.y  = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black"),
          legend.text  = element_text(size = 14, colour = "black"),
          legend.title = element_text(size = 16, colour = "black"))
}
ggsave(file.path(out_dir, "Fig5E_OMG_composition.pdf"),
       plot_composition(c("HAP", "HIF1AKO"), c("WT", "Hif1a KO")), width = 8, height = 7)


## Fig 5F: marker expression in HAP Hif1a KO, on the HAP UMAP
plot_feature <- function(df, gene) {
  ggplot(df, aes(umap_1, umap_2, colour = .data[[gene]])) +
    geom_point(size = 3) +
    scale_colour_gradientn(colours = c("lightgrey", "blue"), name = NULL) +
    ggtitle(gene) + coord_fixed() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
}
genes_5F <- c("Emx2", "Pax6", "Otx1", "Otx2", "Dmbx1", "En1", "En2", "Pax2", "Pax5",   # brain
              "Gata6", "Foxa2", "Sox17")                                               # endoderm
for (g in genes_5F) {
  ggsave(file.path(out_dir, paste0("Fig5F_", g, "_HIF1AKO_FeaturePlot.pdf")),
         plot_feature(ko, g), width = 10, height = 7)
}


## Source Data
xlsx_file <- file.path(sd_dir, "Source_data_Figure_5.xlsx")
if (file.exists(xlsx_file)) file.copy(xlsx_file, sub("\\.xlsx$", "_backup.xlsx", xlsx_file), overwrite = TRUE)
wb <- if (file.exists(xlsx_file)) loadWorkbook(xlsx_file) else createWorkbook()

put_sheet <- function(wb, name, df) {      
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df)
}
put_sheet(wb, "Figure 5E HAP WT UMAP",   hap[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
put_sheet(wb, "Figure 5E Hif1a KO on HAP", ko[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
put_sheet(wb, "Figure 5E bar",           comp_table(c("HAP", "HIF1AKO")))
put_sheet(wb, "Figure 5F",               ko[, c("cell_id", "umap_1", "umap_2", "cell_state", genes_5F)])
saveWorkbook(wb, xlsx_file, overwrite = TRUE)


