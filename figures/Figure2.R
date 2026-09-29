library(ggplot2)
library(dplyr)
library(openxlsx)
source("scripts/OMG_colors.r")

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hap <- read.delim(file.path(sd_dir, "HAP_UMAP_expression.tsv.gz"), check.names = FALSE)
hap$cell_type_num <- factor(unname(celltype_order[hap$cell_state]),
                            levels = unname(celltype_order))

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Fig 2F: HAP UMAP, coloured and numbered by OMG cell state
centroids <- hap %>% group_by(cell_type_num, state_id) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p <- ggplot(hap, aes(umap_1, umap_2, colour = cell_type_num)) +
  geom_point(size = 3) +
  scale_colour_manual(name = "", values = cell_type_colored_numbered) +
  geom_text(data = centroids, aes(umap_1, umap_2, label = state_id),
            size = 6, fontface = "bold", colour = "black") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  coord_fixed() + theme_umap
ggsave(file.path(out_dir, "Fig2F_HAP_UMAP.pdf"), p, width = 10, height = 7)


## Fig 2G: marker expression on the HAP UMAP 
##   grey -> blue = log-normalised expression, scale per gene (FeaturePlot default)
plot_feature <- function(df, gene) {
  ggplot(df, aes(umap_1, umap_2, colour = .data[[gene]])) +
    geom_point(size = 3) +
    scale_colour_gradientn(colours = c("lightgrey", "blue"), name = NULL) +
    ggtitle(gene) +
    coord_fixed() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
}

genes_2G <- c("Emx2", "Pax6", "Otx2", "Otx1", "Dmbx1", "Wnt1",      # forebrain / midbrain
              "Egr2", "Gbx2",                                        # hindbrain
              "En1", "En2", "Pax2", "Pax5",                          # MHB
              "Zic1", "Pax7", "Dbx1", "Shh", "Nkx6-1", "Olig2",      # neural tube
              "Sox10", "Tfap2a")                                     # neural crest

# one PDF per gene
for (g in genes_2G) {
  ggsave(file.path(out_dir, paste0("Fig2G_", g, "_FeaturePlot.pdf")),
         plot_feature(hap, g), width = 10, height = 7)
}


## Source Data: one sheet per panel
wb <- createWorkbook()
addWorksheet(wb, "Fig2F"); writeData(wb, "Fig2F", hap[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
addWorksheet(wb, "Fig2G"); writeData(wb, "Fig2G", hap[, c("cell_id", "umap_1", "umap_2", genes_2G)])
saveWorkbook(wb, file.path(sd_dir, "Source_Data_Fig2.xlsx"), overwrite = TRUE)