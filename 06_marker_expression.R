
# Per-condition marker expression and cross-condition UMAP projections
#
# Input : data/scRNAseq/Asmb_OMG.rds (all four conditions, OMG-annotated)
#         data/scRNAseq/HAP_OMG.rds  (HAP reference embedding, from 04)
#
# Steps : 1. Marker FeaturePlots per condition (HAP, HAP+XAV, HIF1A-KO, NAP),
#            each re-normalized / scaled / PCA / UMAP in isolation
#         2. Per-condition UMAPs coloured by cell state and by predicted day
#         3. Project NAP and HIF1A-KO onto the HAP reference UMAP
#            (FindTransferAnchors + MapQuery) and plot over the HAP background,
#            with projected marker FeaturePlots
#         4. HAP vs HAP+XAV on a shared coordinate space + common colour scale
#         5. Joint HAP / HAP+XAV UMAP (merge, PCA and UMAP) coloured by
#            condition and by cell state, plus joint-UMAP marker FeaturePlots
#
# Output: FeaturePlots in images/scRNAseq/Marker_Expression/,
#         condition / projection UMAPs in images/scRNAseq/OMG/

##############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(ggrastr)
library(tibble)
library(tidyr)
set.seed(42)
source("scripts/OMG_colors.r")

Asmb <- readRDS("data/scRNAseq/Asmb_OMG.rds")
HAP <- readRDS("data/scRNAseq/HAP_OMG.rds")

# selected markers
markers <- c("Foxg1", "Otx2", "Lhx2", "Dbx1", "Wnt1", "En1", "En2", "Pax2", "Pax5",
             "Gbx2", "Krox20", "Egr2", "Mafb", "Hoxa1", "Hoxb1", "Hoxb2",
             "Hoxb3", "Hoxb4", "Hoxd4", "Hoxa5",
             "Hoxc6", "Hoxb6", "Hoxb8", "Hoxc8", "Hoxa9",
             "Hoxa10","Hoxb10","Hoxa11","Hoxb11","Hoxa12","Hoxb12","Hoxa13","Hoxb13",
             "Cdx2","Cdx4",
             "T","Sox2", "Foxa2", "Shh", 
             "Foxd3","Tfap2a","Snai2","Pax3","Sox9","Sox10", "Lmx1a","Lmx1b","Msx1","Msx2")
markers_2 <- c(
  "Foxg1", "Six3", "Lhx2", "Otx1", "Otx2", "Emx2", "Pax6", "Dmbx1", "Fgf8", "Wnt1",
  "En1", "En2", "Pax2", "Pax5", "Pax8", "Egr2", "Mafb", "Gbx2", "Hoxb1", "Pax3",
  "Rfx4", "Zic1", "Olig2", "Dbx1", "Hes5", "Dbx2", "Nkx6-1", "T", "Sox2", "Hes7"
)
markers_endo_meso <- c("Gata6", "Gata4", "Hhex", "Sox17", "Pax9", "Pax1", "Sim1", "Pax7", "Pax3", "EphrinB2", "Tbx18", "EphA4", "Uncx", "Msgn1", "Mesp2", "Meox1", "Foxc2", "Foxc1", "Tbx6")
markers <- c(markers, markers_2, markers_endo_meso)

# check first if all markers are present in the data
markers_present <- markers[markers %in% rownames(HAP)]
# print missing markers
markers_missing <- markers[!markers %in% rownames(HAP)]
cat("Markers present in the data:", markers_present, "\n")
cat("Markers missing in the data:", markers_missing, "\n")

# for each marker, plot feature plot in the umap
Idents(HAP) <- HAP$predicted.celltype_updated
for (marker in markers_present) {
  p <- FeaturePlot(HAP, features = marker, label = TRUE, pt.size = 3, label.size = 2) + 
    ggtitle(marker) + 
    theme(plot.title = element_text(hjust = 0.5)) + # remove axes and ticks
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank()) + coord_fixed(ratio = 1)
    pdf(paste0("images/scRNAseq/Marker_Expression/", marker, "_HAP_FeaturePlot.pdf"), width = 10)
  print(p)
  dev.off()
}

HAPX <- subset(Asmb, subset = condition == "Hypo_XAV")
HIF1 <- subset(Asmb, subset = condition == "HIF1AKO")
NAP <- subset(Asmb, subset = condition == "Normoxic")

HAPX <- NormalizeData(HAPX)
HIF1 <- NormalizeData(HIF1)
NAP <- NormalizeData(NAP)

all.genes <- rownames(HAPX)
HAPX <- ScaleData(HAPX, features=all.genes)
HAPX <- RunPCA(HAPX, npcs=50, features = VariableFeatures(HAPX))
HAPX <- RunUMAP(HAPX, dims = 1:30, return.model=TRUE)

Idents(HAPX) <- HAPX$predicted.celltype_updated
for (marker in markers_present) {
  p <- FeaturePlot(HAPX, features = marker, label = TRUE, pt.size = 3, label.size = 2) + 
    ggtitle(marker) + 
    theme(plot.title = element_text(hjust = 0.5))+ # remove axes and ticks
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank()) + coord_fixed(ratio = 1)
    pdf(paste0("images/scRNAseq/Marker_Expression/", marker, "_HAPX_FeaturePlot.pdf"), width = 10)
  print(p)
  dev.off()
}

all.genes <- rownames(HIF1)
HIF1 <- ScaleData(HIF1, features=all.genes)
HIF1 <- RunPCA(HIF1, npcs=50, features = VariableFeatures(HIF1))
HIF1 <- RunUMAP(HIF1, dims = 1:30, return.model=TRUE)
Idents(HIF1) <- HIF1$predicted.celltype_updated

for (marker in markers_present) {
  p <- FeaturePlot(HIF1, features = marker, label = TRUE, pt.size = 3, label.size = 2) + 
    ggtitle(marker) + 
    theme(plot.title = element_text(hjust = 0.5))+ # remove axes and ticks
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank()) + coord_fixed(ratio = 1)
    pdf(paste0("images/scRNAseq/Marker_Expression/", marker, "_HIF1_FeaturePlot.pdf"), width = 10)
  print(p)
  dev.off()
}

all.genes <- rownames(NAP)
NAP <- ScaleData(NAP, features=all.genes)
NAP <- RunPCA(NAP, npcs=50, features = VariableFeatures(NAP))
NAP <- RunUMAP(NAP, dims = 1:30, return.model=TRUE)
Idents(NAP) <- NAP$predicted.celltype_updated 
for (marker in markers_present) {
  p <- FeaturePlot(NAP, features = marker, label = TRUE, pt.size = 3, label.size = 2) + 
    ggtitle(marker) + 
    theme(plot.title = element_text(hjust = 0.5))+ # remove axes and ticks
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank()) + coord_fixed(ratio = 1)
    pdf(paste0("images/scRNAseq/Marker_Expression/", marker, "_NAP_FeaturePlot.pdf"), width = 10)
  print(p)
  dev.off()
}

# for each object plot umap colored by cell type
umap_data_HAP <- as.data.frame(HAP[["umap"]]@cell.embeddings)
umap_data_HAPX <- as.data.frame(HAPX[["umap"]]@cell.embeddings)
umap_data_HIF1 <- as.data.frame(HIF1[["umap"]]@cell.embeddings)
umap_data_NAP <- as.data.frame(NAP[["umap"]]@cell.embeddings)

umap_data_HAP$condition <- HAP@meta.data[row.names(umap_data_HAP),"condition"]
umap_data_HAPX$condition <- HAPX@meta.data[row.names(umap_data_HAPX),"condition"]
umap_data_HIF1$condition <- HIF1@meta.data[row.names(umap_data_HIF1),"condition"]
umap_data_NAP$condition <- NAP@meta.data[row.names(umap_data_NAP),"condition"]

umap_data_HAP$cell_state_annot <- HAP@meta.data[row.names(umap_data_HAP),"predicted.celltype_updated"]
umap_data_HAPX$cell_state_annot <- HAPX@meta.data[row.names(umap_data_HAPX),"predicted.celltype_updated"]
umap_data_HIF1$cell_state_annot <- HIF1@meta.data[row.names(umap_data_HIF1),"predicted.celltype_updated"]
umap_data_NAP$cell_state_annot <- NAP@meta.data[row.names(umap_data_NAP),"predicted.celltype_updated"]

umap_data_HAP$day <- HAP@meta.data[row.names(umap_data_HAP),"predicted.day"]
umap_data_HAPX$day <- HAPX@meta.data[row.names(umap_data_HAPX),"predicted.day"]
umap_data_HIF1$day <- HIF1@meta.data[row.names(umap_data_HIF1),"predicted.day"]
umap_data_NAP$day <- NAP@meta.data[row.names(umap_data_NAP),"predicted.day"]


umap_combined <- bind_rows(
  umap_data_HAP,
  umap_data_HAPX,
  umap_data_HIF1,
  umap_data_NAP
)
umap_combined$condition <- factor(umap_combined$condition)

unique(umap_combined$cell_state_annot)

umap_combined$cell_state_annot <- factor(
  umap_combined$cell_state_annot,
  levels = names(cell_type_colored)
)
unique(umap_combined$cell_state_annot)
umap_combined$cell_type_num <- unname(
  celltype_order[ as.character(umap_combined$cell_state_annot) ]
)
umap_combined$cell_type_num <- factor(
  umap_combined$cell_type_num,
  levels = unname(celltype_order)   
)
head(umap_combined)
umap_combined$num <- str_extract(umap_combined$cell_type_num, "\\d+")
# Define unique conditions
conditions <- unique(umap_combined$condition)


# Loop through each condition and save a separate plot
for (cond in conditions) {
  
  # Filter data for the current condition
  umap_subset <- subset(umap_combined, condition == cond)
  
  centroids <- umap_subset %>%
    group_by(cell_type_num, num) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

  single_plot <- ggplot() +
    geom_point(data = umap_subset, aes(x = umap_1, y = umap_2, color = cell_type_num), 
               size = 3) +
    scale_color_manual(name = "", values = cell_type_colored_numbered) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # add num for each cell type 
    geom_text(data = centroids, aes(x = umap_1, y = umap_2, label = num),
              size = 6, fontface = "bold", color = "black") +
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
  
  ggsave(paste0("images/scRNAseq/OMG/UMAP_", cond, "_celltype.pdf"),
         single_plot, width = 10, height=7)

  print(cond)
  print(dim(umap_subset))
}

# Loop through each condition and save a separate plot
for (cond in conditions) {
  
  # Filter for current condition
  umap_subset <- subset(umap_combined, condition == cond)
  
  centroids <- umap_subset %>%
    group_by(cell_type_num, num) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

  single_plot <- ggplot() +
    geom_point(data = umap_subset, aes(x = umap_1, y = umap_2, color = day), 
               size = 3) +
    scale_color_manual(name = "", values = day_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # add num for each cell type 
    geom_text(data = centroids, aes(x = umap_1, y = umap_2, label = num),
              size = 6, fontface = "bold", color = "black") +
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
  
  ggsave(paste0("images/scRNAseq/OMG/UMAP_", cond, "_day.pdf"),
         single_plot, width = 10, height = 7)

  print(cond)
  print(dim(umap_subset))
}

#############################################################################################################################################
# UMAP projection of NAP and the KO (HIF1) condition on the HAP umap

# Project each query onto the HAP reference UMAP (coordinates only)
project_onto_hap <- function(query_obj, query_name) {
  anchors <- FindTransferAnchors(
    reference = HAP,
    query = query_obj,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:30
  )
  mapped <- MapQuery(
    anchorset = anchors,
    query = query_obj,
    reference = HAP,
    reference.reduction = "pca",
    reduction.model = "umap"
  )
  proj <- as.data.frame(mapped[["ref.umap"]]@cell.embeddings)
  colnames(proj)[1:2] <- c("umap_1", "umap_2")
  proj$condition <- query_name
  proj$cell_state_annot <- query_obj@meta.data[row.names(proj), "predicted.celltype_updated"]
  proj
}

umap_NAP_proj  <- project_onto_hap(NAP,  "Normoxic")
umap_HIF1_proj <- project_onto_hap(HIF1, "HIF1AKO")

# Plot projected dataset over the HAP background
for (proj_df in list(umap_NAP_proj, umap_HIF1_proj)) {
  cond <- unique(proj_df$condition)

  proj_df$cell_type_num <- unname(celltype_order[as.character(proj_df$cell_state_annot)])
  proj_df$cell_type_num <- factor(proj_df$cell_type_num, levels = unname(celltype_order))
  proj_df$num <- str_extract(proj_df$cell_type_num, "\\d+")

  centroids <- proj_df %>%
    filter(!is.na(cell_type_num)) %>%
    group_by(cell_type_num, num) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

  p <- ggplot() +
    geom_point(data = umap_data_HAP, aes(x = umap_1, y = umap_2),
               color = "grey82", size = 3, alpha = 0.5) +
    geom_point(data = proj_df, aes(x = umap_1, y = umap_2, color = cell_type_num),
               size = 3) +
    scale_color_manual(name = "", values = cell_type_colored_numbered, na.value = "gray80") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_text(data = centroids, aes(x = umap_1, y = umap_2, label = num),
              size = 6, fontface = "bold", color = "black") +
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1)

  ggsave(paste0("images/scRNAseq/OMG/UMAP_", cond, "_on_HAP.pdf"),
         p, width = 10, height = 7)
}

# Feature plot of markers for the projected datasets
map_onto_hap <- function(query_obj) {
  anchors <- FindTransferAnchors(
    reference = HAP,
    query = query_obj,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:30
  )
  MapQuery(
    anchorset = anchors,
    query = query_obj,
    reference = HAP,
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}

mapped_NAP  <- map_onto_hap(NAP)
mapped_HIF1 <- map_onto_hap(HIF1)

proj_entries <- list(
  list(obj = mapped_NAP,  name = "Normoxic"),
  list(obj = mapped_HIF1, name = "HIF1AKO")
)

for (marker in markers_present) {
  for (entry in proj_entries) {
    p <- FeaturePlot(entry$obj, features = marker, reduction = "ref.umap",
                     label = TRUE, pt.size = 3, label.size = 2) +
      ggtitle(marker) +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_blank()) + coord_fixed(ratio = 1)
    ggsave(paste0("images/scRNAseq/OMG/", marker, "_", entry$name, "_projected_FeaturePlot.pdf"),
           p, width = 10, height = 7)
  }
}




#############################################################################################################################################
# Joined UMAP of HAP and HAP+X, colored by condition 

HAP_HAPX <- merge(HAP, y = HAPX, add.cell.ids = c("HAP", "HAPX"))
HAP_HAPX <- JoinLayers(HAP_HAPX)
HAP_HAPX <- NormalizeData(HAP_HAPX)  
# order panels 
HAP_HAPX$condition <- factor(HAP_HAPX$condition, levels = c("Hypoxic", "Hypo_XAV"))
Idents(HAP_HAPX) <- HAP_HAPX$predicted.celltype_updated
HAP_HAPX <- FindVariableFeatures(HAP_HAPX)
HAP_HAPX <- ScaleData(HAP_HAPX, features = rownames(HAP_HAPX))
HAP_HAPX <- RunPCA(HAP_HAPX, npcs = 50)        
HAP_HAPX <- RunUMAP(HAP_HAPX, dims = 1:30, return.model = TRUE)

# plot umap colored by predicted.celltype_updated, but only in one condition at a time (HAP vs HAPX) 
emb <- Embeddings(HAP_HAPX, "umap")
umap_joint <- data.frame(
  umap_1 = emb[, 1],
  umap_2 = emb[, 2],
  condition = HAP_HAPX$condition,
  cell_state_annot = HAP_HAPX$predicted.celltype_updated
)

# shared axis ranges 
xlim <- range(umap_joint$umap_1)
ylim <- range(umap_joint$umap_2)

for (cond in c("Hypoxic", "Hypo_XAV")) {
  umap_subset <- subset(umap_joint, condition == cond)

  umap_subset$cell_type_num <- unname(
    celltype_order[as.character(umap_subset$cell_state_annot)])
  umap_subset$cell_type_num <- factor(umap_subset$cell_type_num,
                                      levels = unname(celltype_order))
  umap_subset$num <- str_extract(umap_subset$cell_type_num, "\\d+")

  centroids <- umap_subset %>%
    filter(!is.na(cell_type_num)) %>%
    group_by(cell_type_num, num) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

  p <- ggplot() +
    # grey background: ALL cells regardless of condition
    geom_point(data = umap_joint,
               aes(x = umap_1, y = umap_2),
               color = "grey82", size = 3, alpha = 0.5) +
    # colored foreground: only this condition
    geom_point(data = umap_subset,
               aes(x = umap_1, y = umap_2, color = cell_type_num),
               size = 3) +
    scale_color_manual(name = "", values = cell_type_colored_numbered,
                       na.value = "gray80") +
    geom_text(data = centroids, aes(x = umap_1, y = umap_2, label = num),
              size = 6, fontface = "bold", color = "black") +
    ggtitle(cond) +
    xlim(xlim) + ylim(ylim) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )

  ggsave(paste0("images/scRNAseq/OMG/",
                "HAP_HAPX_jointUMAP_", cond, "_celltype_numbered.pdf"),
         p, width = 10, height = 7)
}


pdf("images/scRNAseq/OMG/HAP_HAPX_jointUMAP.pdf",
    width = 10, height = 7)
DimPlot(HAP_HAPX, reduction = "umap", group.by = "condition", pt.size = 3) &
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), axis.line = element_blank()) &
  coord_fixed(ratio = 1)
dev.off()

#############################################################################################################################################
# Joint-UMAP marker FeaturePlots, per condition, shared colour range

markers_shared <- markers[markers %in% rownames(HAP_HAPX)]

# shared axis ranges from the full joint umap 
emb  <- Embeddings(HAP_HAPX, "umap")
xlim <- range(emb[, 1]); ylim <- range(emb[, 2])

for (marker in markers_shared) {
  # global max across both conditions 
  vmax <- max(GetAssayData(HAP_HAPX, layer = "data")[marker, ], na.rm = TRUE)

  for (cond in c("Hypoxic", "Hypo_XAV")) {
    obj <- subset(HAP_HAPX, subset = condition == cond)

    p <- FeaturePlot(obj, features = marker, reduction = "umap",   # joint umap
                     pt.size = 3) +
      scale_color_gradientn(colours = c("lightgrey", "blue"),
                            limits = c(0, vmax),
                            oob = scales::squish) +
      ggtitle(paste(marker, cond)) +
      xlim(xlim) + ylim(ylim) +                 # shared axes, both panels match
      coord_fixed(ratio = 1) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title = element_blank(), axis.line = element_blank())

    ggsave(paste0("images/scRNAseq/Marker_Expression/",
                  marker, "_", cond, "_jointUMAP_FeaturePlot.pdf"),
           p, width = 10, height = 7)
  }
}
