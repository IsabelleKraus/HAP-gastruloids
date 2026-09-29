library(ggplot2)
library(dplyr)
library(openxlsx)
source("scripts/OMG_colors.r")

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## Fig 3C (top): Qiu 2024 reference-predicted somite count per HAP cell, mean +/- SD
som <- read.delim(file.path(sd_dir, "HAP_predicted_somite_count.tsv.gz"))

df_sum <- data.frame(group = "HAP",
                     mean  = mean(som$predicted_somite_count, na.rm = TRUE),
                     sd    = sd(som$predicted_somite_count,   na.rm = TRUE))

set.seed(42)   # fixed jitter positions
p <- ggplot() +
  geom_col(data = df_sum, aes(x = group, y = mean), fill = "#963752", width = 0.6) +
  geom_errorbar(data = df_sum, aes(x = group, ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_jitter(data = som, aes(x = "HAP", y = predicted_somite_count),
              width = 0.15, size = 0.5, alpha = 0.3, colour = "black") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(y = "Predicted somite count (HAPs)") +
  theme_classic() +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 16))
ggsave(file.path(out_dir, "Fig3C_predicted_somite_count.pdf"), p, width = 2.5, height = 5)


## Fig 3D: HAP UMAP + (paraxial) mesoderm markers 
hap <- read.delim(file.path(sd_dir, "HAP_UMAP_expression.tsv.gz"), check.names = FALSE)
hap$cell_type_num <- factor(unname(celltype_order[hap$cell_state]),
                            levels = unname(celltype_order))

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

centroids <- hap %>% group_by(cell_type_num, state_id) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
p <- ggplot(hap, aes(umap_1, umap_2, colour = cell_type_num)) +
  geom_point(size = 3) +
  scale_colour_manual(name = "", values = cell_type_colored_numbered) +
  geom_text(data = centroids, aes(umap_1, umap_2, label = state_id),
            size = 6, fontface = "bold", colour = "black") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  coord_fixed() + theme_umap
ggsave(file.path(out_dir, "Fig3D_HAP_UMAP.pdf"), p, width = 10, height = 7)

plot_feature <- function(df, gene) {
  ggplot(df, aes(umap_1, umap_2, colour = .data[[gene]])) +
    geom_point(size = 3) +
    scale_colour_gradientn(colours = c("lightgrey", "blue"), name = NULL) +
    ggtitle(gene) + coord_fixed() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
}
genes_3D <- c("Foxc1", "Foxc2", "Hes7", "Meox1", "Mesp2", "Pax1", "Pax3", "Pax9",
              "T", "Tbx6", "Tbx18", "Uncx")
for (g in genes_3D) {
  ggsave(file.path(out_dir, paste0("Fig3D_", g, "_FeaturePlot.pdf")),
         plot_feature(hap, g), width = 10, height = 7)
}


## Source Data: 
xlsx_file <- file.path(sd_dir, "Source data Figure 3.xlsx")
file.copy(xlsx_file, sub("\\.xlsx$", "_backup.xlsx", xlsx_file), overwrite = TRUE)
wb <- loadWorkbook(xlsx_file)

# 3C top: per-cell values (A-C) + the plotted summary (E-F)
writeData(wb, "Figure 3C top", som, startCol = 1, startRow = 1)
writeData(wb, "Figure 3C top",
          data.frame(statistic = c("mean", "SD", "n cells"),
                     value     = c(df_sum$mean, df_sum$sd, nrow(som))),
          startCol = 5, startRow = 1)

# 3D: UMAP coordinates, cell state and marker expression per cell
if (!"Figure 3D" %in% names(wb)) addWorksheet(wb, "Figure 3D")
writeData(wb, "Figure 3D", hap[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id", genes_3D)])

saveWorkbook(wb, xlsx_file, overwrite = TRUE)
