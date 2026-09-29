library(ggplot2)
library(dplyr)
library(openxlsx)
source("scripts/OMG_colors.r")    # cell_type_colored_numbered, celltype_order

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hap <- read.delim(file.path(sd_dir, "HAP_UMAP_expression.tsv.gz"), check.names = FALSE)

# numbered labels, e.g. "Midbrain (3)"; factor levels fix the legend order
hap$cell_type_num <- factor(unname(celltype_order[hap$cell_state]),
                            levels = unname(celltype_order))

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## ED Fig 5E (left): HAP UMAP, coloured and numbered by OMG cell state
centroids <- hap %>% group_by(cell_type_num, state_id) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
p <- ggplot(hap, aes(umap_1, umap_2, colour = cell_type_num)) +
  geom_point(size = 3) +
  scale_colour_manual(name = "", values = cell_type_colored_numbered) +
  geom_text(data = centroids, aes(umap_1, umap_2, label = state_id),
            size = 6, fontface = "bold", colour = "black") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  coord_fixed() + theme_umap
ggsave(file.path(out_dir, "ED5E_HAP_UMAP.pdf"), p, width = 10, height = 7)


## ED Fig 5E (right): endoderm marker expression on the HAP UMAP
##   expression, colour scale per gene (FeaturePlot default)
plot_feature <- function(df, gene) {
  ggplot(df, aes(umap_1, umap_2, colour = .data[[gene]])) +
    geom_point(size = 3) +
    scale_colour_gradientn(colours = c("lightgrey", "blue"), name = NULL) +
    ggtitle(gene) + coord_fixed() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
}
genes_ED5E <- c("Cdx2", "Foxa2", "Gata4", "Gata6", "Hhex", "Sox17")
for (g in genes_ED5E) {
  ggsave(file.path(out_dir, paste0("ED5E_", g, "_HAP_FeaturePlot.pdf")),
         plot_feature(hap, g), width = 10, height = 7)
}


## Source Data
xlsx_file <- file.path(sd_dir, "Source_data_Extended_Data_Figure_5.xlsx")
if (file.exists(xlsx_file)) file.copy(xlsx_file, sub("\\.xlsx$", "_backup.xlsx", xlsx_file), overwrite = TRUE)
wb <- if (file.exists(xlsx_file)) loadWorkbook(xlsx_file) else createWorkbook()

put_sheet <- function(wb, name, df) {         
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df)
}
put_sheet(wb, "ED Fig 5E HAP UMAP", hap[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
put_sheet(wb, "ED Fig 5E markers",  hap[, c("cell_id", "umap_1", "umap_2", "cell_state", genes_ED5E)])
saveWorkbook(wb, xlsx_file, overwrite = TRUE)