library(ggplot2)
library(dplyr)
library(ggrastr)
library(openxlsx)
source("scripts/OMG_colors.r")    # cell_type_colored(_numbered), celltype_order

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref   <- read.delim(file.path(sd_dir, "OMG_reference_UMAP.tsv.gz"))
query <- read.delim(file.path(sd_dir, "OMG_query_UMAP.tsv.gz"))

# numbered labels, e.g. "Midbrain (3)"; factor levels fix the legend order
add_num <- function(df) {
  df$cell_type_num <- factor(unname(celltype_order[df$cell_state]),
                             levels = unname(celltype_order))
  df
}
ref   <- add_num(ref)
query <- add_num(query)

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## ED Fig 9E (left): HAP and HAP + XAV projected onto the OMG reference

ref_rep <- ref[rep(seq_len(nrow(ref)), 5), ]

plot_projection <- function(cond) {
  sub <- subset(query, condition == cond)
  centroids <- sub %>% group_by(cell_type_num, state_id) %>%
    summarise(refumap_1 = median(refumap_1), refumap_2 = median(refumap_2), .groups = "drop")
  ggplot() +
    rasterise(geom_point(data = ref_rep, aes(umap_1, umap_2),
                         size = 3, colour = "grey82", alpha = 0.5), dpi = 300) +
    geom_point(data = sub, aes(refumap_1, refumap_2, colour = cell_type_num), size = 3) +
    scale_colour_manual(name = "", values = cell_type_colored_numbered) +
    geom_text(data = centroids, aes(refumap_1, refumap_2, label = state_id),
              size = 6, fontface = "bold", colour = "black") +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    coord_fixed() + theme_umap
}
ggsave(file.path(out_dir, "ED9E_UMAP_HAP_projected_OMG.pdf"),  plot_projection("HAP"),  width = 10, height = 7)
ggsave(file.path(out_dir, "ED9E_UMAP_HAPX_projected_OMG.pdf"), plot_projection("HAPX"), width = 10, height = 7)

## inset: OMG reference coloured by cell state, no labels / legend (small version of Fig. 1G)
p <- ggplot(ref, aes(umap_1, umap_2, colour = cell_type_num)) +
  rasterise(geom_point(size = 0.5), dpi = 300) +
  scale_colour_manual(values = cell_type_colored_numbered) +
  coord_fixed() + theme_void() + theme(legend.position = "none")
ggsave(file.path(out_dir, "ED9E_inset_OMG_reference.pdf"), p, width = 3, height = 3)


## ED Fig 9E (right): composition bar HAP vs HAP + XAV
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
ggsave(file.path(out_dir, "ED9E_OMG_composition.pdf"),
       plot_composition(c("HAP", "HAPX"), c("HAP", "HAP + XAV")), width = 8, height = 7)


## Source Data: add the scRNA-seq sheets to the Extended Data Figure 9 workbook
##   (opens the existing file if already started; other sheets stay untouched)
xlsx_file <- file.path(sd_dir, "Source_data_Extended_Data_Figure_9.xlsx")
if (file.exists(xlsx_file)) file.copy(xlsx_file, sub("\\.xlsx$", "_backup.xlsx", xlsx_file), overwrite = TRUE)
wb <- if (file.exists(xlsx_file)) loadWorkbook(xlsx_file) else createWorkbook()

put_sheet <- function(wb, name, df) {         # (re)write one sheet
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df)
}
put_sheet(wb, "ED Fig 9E OMG reference", ref[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
put_sheet(wb, "ED Fig 9E HAP HAPX on OMG",
          subset(query, condition %in% c("HAP", "HAPX"),
                 select = c(cell_id, condition, refumap_1, refumap_2, cell_state, state_id)))
put_sheet(wb, "ED Fig 9E bar", comp_table(c("HAP", "HAPX")))
saveWorkbook(wb, xlsx_file, overwrite = TRUE)