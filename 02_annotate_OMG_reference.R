
# Refine the OMG reference (E8.0-E9.75)
#
# Input : data/OMG/OMG_E8_E10_merge.rds  (from 01_build_OMG_reference.R)
#
# Steps : 1. Drop E10.0 
#            (data/OMG/OMG_E8_E9.75_merge.rds)
#         2. Subcluster "Facial mesenchyme" (FindClusters res 0.1) into
#            FM_1/FM_2/FM_3 and write these back as celltype_updated
#         3. Rename four states (Diencephalon -> Posterior Forebrain /
#            Diencephalon; Dermomyotome -> Somites / Dermomyotome; Eye field ->
#            Anterior Forebrain; Sclerotome -> Somites / Sclerotome)
#         4. save the annotated reference with the celltype_updated column added
#
# Output: data/OMG/OMG_E8_E9.75_merge.rds (annotated; input to 04)
#         data/markers/OMG_celltype_markers.csv, FM_subclustered_markers.csv
#         figures in images/scRNAseq/OMG/
 
###############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(ggrastr)
source("scripts/OMG_colors.r")
source("data/OMG/JAX_color_code.R")
set.seed(42)

OMG <- readRDS('data/OMG/OMG_E8_E10_merge.rds')
OMG <- subset(OMG, day!="E10.0")

# Find markers for each cell state 
Idents(OMG) <- OMG$celltype_update
markers <- FindAllMarkers(OMG, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
# filter padj < 0.05
markers <- markers %>% filter(p_val_adj < 0.05 & avg_log2FC > 1)
# get top 50 markers per cluster
markers <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.csv(markers, file = 'data/markers/OMG_celltype_markers.csv')
markers <- read.csv('data/markers/OMG_celltype_markers.csv')

# plot heatmap showing top50 marker expression as zscore grouped by cell type (color at the row instead of rownames) and columns order cells by day 
marker_genes <- unique(markers$gene)
unique(OMG$celltype_update)
rownames(OMG[["RNA"]]@scale.data)

# plot heatmap but leave out the following idents:
extraembryonic_states <- c(
  'Extraembryonic visceral endoderm',
  'Amniotic ectoderm'
)
markers <- subset(markers, cluster %in% unique(OMG$celltype_update))
marker_genes <- unique(markers$gene)
present_genes <- marker_genes[marker_genes %in% rownames(OMG[["RNA"]]@scale.data)]
length(present_genes) # Should be >0

pdf('images/OMG_markers/Heatmap_OMG_cellstate_markers.pdf', width = 10, height = 18)
DoHeatmap(OMG, features = present_genes, raster = FALSE)
dev.off()


# plot dotplot for each marker gene per cell type and day 
Idents(OMG) <- OMG$celltype_update
pdf('images/OMG_markers/Dotplot_OMG_cellstate_markers.pdf', width = 20, height = 10)
DotPlot(OMG, features = marker_genes, group.by = 'day', split.by = 'ident') + RotatedAxis()
dev.off()

# complex heatmap plot scaled expression values and ad column annotation for day and cell type
OMG <- ScaleData(OMG, features = unique(marker_genes))
head(OMG[["RNA"]]@scale.data[marker_genes, 1:5])
# extract scaled expression data for marker genes
expr_matrix <- OMG[["RNA"]]@scale.data[marker_genes, ]
# create a data frame for column annotations
col_anno <- data.frame(
  Day = OMG$day,
  CellType = OMG$celltype_update
)
rownames(col_anno) <- colnames(expr_matrix)
# define colors for annotations
unique(col_anno$Day)    
day_colors <- c("E8.0-E8.5" = "#1f77b4", "E8.75" = "#ff7f0e", "E9.0" = "#2ca02c", "E9.25" = "#d62728", "E9.5" = "#9467bd", "E9.75" = "#8c564b", "E10.0" = "#e377c2")

# set colors for cell types autonmatically
cell_types <- unique(col_anno$CellType)
length(cell_types)

# generate distinct colors for each cell type
set.seed(42) 
cell_type_colors <- setNames(rainbow(length(cell_types)), cell_types)

col_colors <- list(
    CellType = cell_type_colors,
    Day = day_colors
)

#first group columns by cell type then by day
column_order <- order(col_anno$CellType, col_anno$Day)
col_anno <- col_anno[column_order, ]

ha <- HeatmapAnnotation(
  df = col_anno,
  col = col_colors,
  which = "column"
)

# order rows genes in expr_matrix by cell type (cluster) in marker, same cell state order as in columns
cluster_order <- unique(col_anno$CellType)
ordered_markers <- markers %>%
  filter(gene %in% rownames(expr_matrix)) %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  arrange(cluster, desc(avg_log2FC))
gene_order <- unique(ordered_markers$gene)

# colors for z score
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf('images/OMG_markers/Heatmap_OMG_cellstate_markers_noraster.pdf', width = 15, height = 40)
Heatmap(expr_matrix[gene_order , column_order], 
        name = "Expression",
        top_annotation = ha,
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 1),
        column_title = "OMG Cell State Marker Genes",
        heatmap_legend_param = list(title = "Scaled Expression"),
        use_raster = FALSE,
        col = col_fun
)
dev.off()

# save as png 
png('images/OMG_markers/Heatmap_OMG_cellstate_markers.png', width = 15, height = 40, units = 'in', res = 300)
Heatmap(expr_matrix[gene_order , column_order], 
        name = "Expression",
        top_annotation = ha, 
        show_row_names = TRUE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 2),
        column_title = "OMG Cell State Marker Genes",
        heatmap_legend_param = list(title = "Scaled Expression"),
        use_raster = FALSE,
        col = col_fun
)
dev.off()

# Get the universe of gene names from the active assay
active_assay <- DefaultAssay(OMG)
gene_universe <- rownames(GetAssayData(OMG, assay = active_assay, slot = "data"))
if (is.null(gene_universe) || length(gene_universe) == 0) {
  gene_universe <- rownames(GetAssayData(OMG, assay = active_assay, slot = "counts"))
}

# Ensure output dir exists
out_dir <- "images/OMG_markers"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (cell_type in cluster_order) {
  cat("Processing cell type:", cell_type, "\n")

  # subset markers for this cell type
  cell_type_markers <- markers %>% filter(cluster == cell_type)
  cell_type_genes    <- unique(cell_type_markers$gene)
  present_genes      <- intersect(cell_type_genes, gene_universe)

  if (length(present_genes) == 0) {
    cat("No marker genes found (present in object) for:", cell_type, "\n")
    next
  }

  # pick cells of this type 
  cols_subset <- col_anno[col_anno$CellType == cell_type, , drop = FALSE]
  if (nrow(cols_subset) == 0) {
    cat("No cells in col_anno for cell type:", cell_type, "\n")
    next
  }

  # intersect with actual cells
  valid_cells <- intersect(rownames(cols_subset), colnames(OMG))
  if (length(valid_cells) == 0) {
    cat("No overlapping cells between col_anno and OMG for:", cell_type, "\n")
    next
  }
  cols_subset <- cols_subset[valid_cells, , drop = FALSE]

# Fetch expression for the genes across ALL cells first, then subset
 expr_df <- FetchData(OMG, vars = present_genes)         
expr_df <- expr_df[valid_cells, , drop = FALSE]

# Make matrix genes × cells 
expr_mat <- t(as.matrix(expr_df))                     

# Z-score per gene (row-wise)
expr_z <- t(scale(t(expr_mat)))                          

# Remove genes that are all NA
expr_z <- expr_z[rowSums(!is.na(expr_z)) > 0, , drop = FALSE]
if (nrow(expr_z) == 0) { cat("All genes NA after scaling for:", cell_type, "\n"); next }

# Order columns (cells) by Day 
ord <- order(cols_subset$Day)
cols_subset <- cols_subset[ord, , drop = FALSE]
expr_z      <- expr_z[, rownames(cols_subset), drop = FALSE]  # <- now valid (columns are cells)

# Order genes by avg_log2FC (descending), keeping only present genes
ordered_markers <- cell_type_markers |>
  dplyr::filter(gene %in% rownames(expr_z)) |>
  dplyr::arrange(desc(avg_log2FC))
gene_order <- unique(ordered_markers$gene)
if (length(gene_order) == 0) { cat("No ordered markers for:", cell_type, "\n"); next }

rng <- range(expr_z[gene_order, , drop = FALSE], na.rm = TRUE)
col_fun_dynamic <- circlize::colorRamp2(c(rng[1], 0, rng[2]), c("blue", "white", "red"))

  ha <- HeatmapAnnotation(
    df = data.frame(Day = cols_subset$Day, row.names = rownames(cols_subset)),
    col = list(Day = day_colors),
    which = "column"
  )

  cat("gene_order length:", length(gene_order), "\n")
  cat("expr_z dim:", dim(expr_z), "\n")
  cat("Valid cells:", ncol(expr_z), "\n")
  cat("Value range:", paste(rng, collapse = " "), "\n")

  safe_cluster_name <- gsub("[/\\:]", "_", cell_type)
  pdf(file.path(out_dir, paste0("Heatmap_OMG_", safe_cluster_name, "_markers.pdf")),
      width = 15, height = 20)

ht <- ComplexHeatmap::Heatmap(
  expr_z[gene_order, , drop = FALSE],   # genes × cells
  name = "Expression (z)",
  top_annotation = ha,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 6),
  column_title = paste("Marker Genes for", cell_type),
  heatmap_legend_param = list(title = "Scaled Expression"),
  use_raster = FALSE,
  col = col_fun_dynamic
)

ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

  dev.off()

# also plot log-normalized expression values (not scaled)
  expr_subset <- FetchData(OMG, vars = present_genes)
  expr_subset <- expr_subset[valid_cells, , drop = FALSE]
  expr_subset <- t(as.matrix(expr_subset))  # genes × cells
  expr_subset <- expr_subset[rownames(expr_z), , drop = FALSE]  # keep same genes as in z-scored
  expr_subset <- expr_subset[, rownames(cols_subset), drop = FALSE]  # keep same cells as in z-scored
  rng2 <- range(expr_subset[gene_order, , drop = FALSE], na.rm = TRUE)
  col_fun2 <- circlize::colorRamp2(c(rng2[1], 0, rng2[2]), c("blue", "white", "red"))
  pdf(file.path(out_dir, paste0("Heatmap_OMG_", safe_cluster_name, "_markers_lognorm.pdf")),
      width = 15, height = 20)
ht2 <- ComplexHeatmap::Heatmap(
  expr_subset[gene_order, , drop = FALSE],   # genes × cells
  name = "Expression (lognorm)",
  top_annotation = ha,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 6),
  column_title = paste("Marker Genes for", cell_type),
  heatmap_legend_param = list(title = "Log-normalized Expression"),
  use_raster = FALSE,
  col = col_fun2
) 
ComplexHeatmap::draw(ht2, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()   

}


# Process
OMG <- NormalizeData(OMG, verbose = FALSE)
OMG <- FindVariableFeatures(OMG, verbose = FALSE)
OMG <- ScaleData(OMG, verbose = FALSE)
OMG <- RunPCA(OMG, npcs = 50, verbose = FALSE)
OMG <- RunUMAP(OMG, dims = 1:50, return.model = T) #, n.components = 3, min.dist = 0.75, 

day_colors <- day_color_plate[c(3:9)]
names(day_colors)[1] <- 'E8.0-E8.5'

pdf('images/OMG_wo10/OMG_E8_E975_merge_days_t.pdf')
DimPlot(OMG, reduction = "umap", group.by='day', raster=T, cols=day_colors)
dev.off()

pdf('images/OMG_wo10/OMG_E8_E975_merge_days_split.pdf', width=56)
DimPlot(OMG, reduction = "umap", group.by='day', split.by='day', raster=T, cols=day_colors)
dev.off()

## UMAP by trajectory
pdf('images/OMG_wo10/OMG_E8_E975_merge_trajectory_split.pdf', height=9, width=14)
DimPlot(OMG, reduction = "umap", group.by='major_trajectory', label=T, raster=T, cols=major_trajectory_color_plate[levels(factor(OMG@meta.data$major_trajectory))])
dev.off()

## Cell types per day
pdf('images/OMG_wo10/OMG_E8_E975_cell_types_per_day.pdf', width=14)
ggplot(OMG@meta.data, aes(x=celltype_update, fill=day)) + geom_bar() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=day_colors)
dev.off()


# Facial Mesenchyme subclusterung 
markers <- read.csv('data/markers/OMG_celltype_markers.csv')
fm <- subset(markers, cluster == "Facial mesenchyme")

pdf("images/OMG_markers/OMG_E8_E975_merge_Facial_mesenchyme_markers_heatmap.pdf", width=30, height=15)
DoHeatmap(OMG, features = fm$gene, group.by = "celltype_update", group.colors=cell_type_colors, raster=F, size=3, angle=90)
dev.off()

# FM subclustering
subset_fm <- subset(OMG, subset = celltype_update %in% c("Facial mesenchyme"))
# process
subset_fm <- NormalizeData(subset_fm, verbose = FALSE)
subset_fm <- FindVariableFeatures(subset_fm, verbose = FALSE)
subset_fm <- ScaleData(subset_fm, verbose = FALSE)
subset_fm <- RunPCA(subset_fm, npcs = 50, verbose = FALSE)
subset_fm <- RunUMAP(subset_fm, dims = 1:50, return.model = T) 
subset_fm <- FindNeighbors(subset_fm, dims = 1:50)
subset_fm <- FindClusters(subset_fm, resolution = 0.1)

# plot distribution of celltype_update per cluster in a stacked barplot
plot_data <- subset_fm@meta.data

# per seurat cluster plot amount of cells per celltype_update
plot_data_summary <- plot_data %>%
  group_by(seurat_clusters, celltype_update) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(seurat_clusters) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  arrange(seurat_clusters, desc(percentage))

# plot 
pdf("images/scRNAseq/OMG/stacked_FM_subclustered_celltype_update_FM.pdf", width=8, height=6)
ggplot(plot_data_summary, aes(x = factor(seurat_clusters), y = count, fill = celltype_update)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors) + 
  labs(
    x = "Subclusters",
    y = "Cells",
    fill = "Cell states"
  ) 
dev.off()

## compute markers for the clusters
Idents(subset_fm) <- subset_fm$seurat_clusters
fm_markers <- FindAllMarkers(subset_fm, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
# filter markers for adjusted p value < 0.05
fm_markers <- subset(fm_markers, p_val_adj < 0.05 & avg_log2FC > 1)
fm_markers$cell_type_updated <- ifelse(fm_markers$cluster == 0, "FM_0", 
                                 ifelse(fm_markers$cluster == 1, "FM_1", "FM_2"))
write.csv(fm_markers, file="tables/scRNAseq/OMG/FM_subclustered_markers.csv")

## Do heatmap of top 10 markers per cluster
top10 <- fm_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("images/scRNAseq/OMG/FM_subclustered_top10_markers_heatmap.pdf", width=10, height=8)
DoHeatmap(subset_fm, features = top10$gene, group.by = "seurat_clusters", raster=F, size=3, angle=90)
dev.off()

## show top 20 markers per cluster 2 4 and 6
top20 <- fm_markers %>% filter(cluster %in% c(2,4,6)) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
pdf("images/scRNAseq/OMG/FM_markers_clusters_heatmap.pdf", width=10, height=8)
DoHeatmap(subset_fm, features = top20$gene, group.by = "seurat_clusters", raster=F, size=3, angle=90)
dev.off() 

# FM_0, FM1, FM2
subset_fm$FM_subcluster <- paste0("FM_", as.numeric(subset_fm$seurat_clusters)-1)

# in OMG 
fm_vec <- setNames(subset_fm@meta.data$FM_subcluster, Cells(subset_fm))
OMG <- AddMetaData(OMG, fm_vec, col.name = "FM_subcluster")
OMG@meta.data$celltype_updated <- ifelse(
  OMG@meta.data$celltype_update == "Facial mesenchyme",
  OMG@meta.data$FM_subcluster,
  OMG@meta.data$celltype_update)


# Plot dimplot of the new annotation
pdf('images/scRNAseq/OMG/UMAP_Filtered_5_OMG_E8_E975_merge_FM_subclustered_celltypes.pdf', width=14)
DimPlot(OMG, reduction = "umap", group.by='celltype_updated', label=T, raster=F, cols=cell_type_colored)
dev.off()

OMG$celltype_updated <- ifelse(OMG$celltype_updated == "Diencephalon", "Posterior Forebrain / Diencephalon", as.character(OMG$celltype_updated))
OMG$celltype_updated <- ifelse(OMG$celltype_updated == "Dermomyotome", "Somites / Dermomyotome", as.character(OMG$celltype_updated))
OMG$celltype_updated <- ifelse(OMG$celltype_updated == "Eye field", "Anterior Forebrain", as.character(OMG$celltype_updated))
OMG$celltype_updated <- ifelse(OMG$celltype_updated == "Sclerotome", "Somites / Sclerotome", as.character(OMG$celltype_updated))

saveRDS(OMG, 'data/OMG/OMG_E8_E9.75_merge.rds')
