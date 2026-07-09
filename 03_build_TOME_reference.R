
# Build the TOME ("Trajectories Of Mouse Embryogenesis") reference (E7.25-E10.5)
# and stage / annotate the HAP and TLS gastruloids against it
#
# Input : data/TOME/seurat_object_E7.25.rds till E10.5.rds (per-stage TOME
#         objects); data/scRNAseq/Asmb.rds (from 04, QC'd gastruloids) and
#         data/TLS/TLS_norm.Robj for the projection step
#
# Steps : 1. Load per-stage TOME objects, merge, map Ensembl IDs to symbols
#         2. Per-batch normalize / HVG (3,000) / scale / PCA, then RPCA
#            integration (FindIntegrationAnchors + IntegrateData); PCA + UMAP
#         3. Subcluster Hindbrain (res 0.01) and Forebrain/midbrain (res 0.05);
#            annotate subclusters from markers -> cell_type_updated
#            -> save data/TOME/TOME_E7.25_E10.5_rpca.rds
#         4. Project HAP + TLS onto TOME (FindTransferAnchors / MapQuery);
#            transfer cell_type and day; plot prediction scores
#         5. Filter by prediction score (<20th percentile) and drop
#            condition/state pairs with <10 cells
#         6. Predicted-stage dot plots,
#            reference and projected UMAPs, composition bar plots, and a
#            score-weighted expected stage per cell type for HAP
#         7. Filtered-reference markers -> TOME_filtered_markers.csv + heatmap
#
# Output: data/TOME/TOME_E7.25_E10.5_rpca.rds (annotated reference)
#         data/scRNAseq/{TOME_filtered,Asmb_projected_filtered,
#           TLS_projected_filtered}_20_TOME_E725_E105.rds
#         tables/scRNAseq/TOME/*, figures in images/scRNAseq/TOME/ and
#         images/new_strategy/
#
# Note  : Depends on data/scRNAseq/Asmb.rds from 04. The HAP staging result
#         from step 6 is what motivates the E8.0-E9.75 OMG restriction in 02,
#         so in practice: 04 (QC) -> 03 -> 02 -> 04 (OMG annotation).

###############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

library(Seurat)
library(future)
library(future.apply)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(readr)
library(ggrastr)
library(openxlsx)
source("scripts/TOME_colors.r")
set.seed(42)

# load data
E7.25 <- readRDS('data/TOME/seurat_object_E7.25.rds')
E7.5 <- readRDS('data/TOME/seurat_object_E7.5.rds')
E7.75 <- readRDS('data/TOME/seurat_object_E7.75.rds')
E8 <- readRDS('data/TOME/seurat_object_E8.rds')
E8.25 <- readRDS('data/TOME/seurat_object_E8.25.rds')
E8.5a <- readRDS('data/TOME/seurat_object_E8.5a.rds')
E8.5b <- readRDS('data/TOME/seurat_object_E8.5b.rds')
E9.5 <- readRDS('data/TOME/seurat_object_E9.5.rds')
E10.5 <- readRDS('data/TOME/seurat_object_E10.5.rds')

########################################################################################
# # integrate
obj.list <- list(E7.25, E7.5, E7.75, E8, E8.25, E8.5a, E8.5b, E9.5, E10.5)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)

names(obj.list) <- c("E7.25","E7.5","E7.75",
                     "E8","E8.25","E8.5a","E8.5b","E9.5","E10.5")
for (i in seq_along(obj.list)) {
  obj.list[[i]]$stage <- names(obj.list)[i]
}
obj <- Reduce(function(x, y) merge(x, y), obj.list)

## rename genes
library(AnnotationDbi)
library(org.Mm.eg.db)

### Get the counts and strip Ensembl version suffixes
assay_name <- DefaultAssay(obj)
counts <- GetAssayData(obj, assay = assay_name, slot = "counts")
ens_ids <- gsub("\\..*$", "", rownames(counts))

### Map Ensembl -> Symbol
symbols <- mapIds(
  org.Mm.eg.db,
  keys = ens_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

### Fallback: if no symbol, keep the (cleaned) Ensembl ID
symbols_clean <- ifelse(is.na(symbols) | symbols == "", ens_ids, symbols)

### Apply the new gene names and collapse duplicates (same symbol from multiple Ensembl IDs)
rownames(counts) <- symbols_clean

### aggregate.Matrix preserves sparsity; group by rownames and sum
library(Matrix)

#### counts: dgCMatrix (rows = features, cols = cells)
grouping <- rownames(counts)              # target names after mapping (e.g., symbols)
uniq <- unique(grouping)
row_index <- match(grouping, uniq)

#### Build a sparse mapping matrix G (n_groups x n_rows), with 1 at (group, row)
G <- sparseMatrix(
  i = row_index,
  j = seq_along(grouping),
  x = 1,
  dims = c(length(uniq), length(grouping)),
  dimnames = list(uniq, NULL)
)

#### Collapse rows: (n_groups x n_rows) %*% (n_rows x n_cells) -> (n_groups x n_cells)
counts_collapse <- G %*% counts
rownames(counts_collapse) <- uniq

### Rebuild the assay to keep everything consistent
new_assay <- CreateAssayObject(counts = counts_collapse)
obj[[assay_name]] <- new_assay
DefaultAssay(obj) <- assay_name


## process
# split it by batch
objs <- SplitObject(obj, split.by = "orig.ident")

# Per-object preprocessing
objs <- lapply(objs, function(x){
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 3000, verbose = FALSE)
  return(x)
})

# Pick shared features across all objects
features <- SelectIntegrationFeatures(objs, nfeatures = 3000)

# Scale & PCA per object using only these features
objs <- lapply(objs, function(x){
  x <- ScaleData(x, features = features, verbose = FALSE,
                 block.size = 1000)   # smaller blocks = lower RAM spikes
  x <- RunPCA(x, features = features, npcs = 50, verbose = FALSE)
  return(x)
})


# Select features and find anchors using reciprocal PCA
anchors  <- FindIntegrationAnchors(
  object.list = objs,
  anchor.features = features,
  reduction = "rpca",      
  dims = 1:30,            
  k.anchor = 5 ,           
  k.filter = 50 
)

integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Downstream on the integrated assay
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:50, return.model = TRUE)

integrated@meta.data$day <- factor(integrated@meta.data$day, levels = c("E6.75", "E7","E7.25","E7.5","E7.75", "E8","E8.25","E8.5a", "E8.5b", "E9.5","E10.5","E11.5"))

## UMAP by time point
pdf('images/new_strategy/TOME_rpca_E7.25_E10.5_days.pdf')
DimPlot(integrated, reduction = "umap", group.by='day', raster=F)
dev.off()

pdf('images/new_strategy/TOME_rpca_E7.25_E10.5_celltype.pdf', width=14)
DimPlot(integrated, reduction = "umap", group.by='cell_type', raster=F)
dev.off()

# save integrated
saveRDS(integrated, 'data/TOME/TOME_E7.25_E10.5_rpca.rds')

# hindbrain subclustering 
hb <- subset(integrated, subset = cell_type %in% c("Hindbrain"))

# subcluster processing
hb <- ScaleData(hb, verbose = FALSE)
hb <- RunPCA(hb, npcs = 50, verbose = FALSE)
hb <- FindNeighbors(hb, dims = 1:50)
hb <- FindClusters(hb, resolution = 0.01)

plot_data <- hb@meta.data

# per seurat cluster plot amount of cells per celltype_update
plot_data_summary <- plot_data %>%
  group_by(seurat_clusters, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(seurat_clusters) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  arrange(seurat_clusters, desc(percentage))

# plot 
pdf("images/scRNAseq/TOME/bar_Hindbrain_subclustered.pdf", width=8, height=6)
ggplot(plot_data_summary, aes(x = factor(seurat_clusters), y = count, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors) + 
  labs(
    x = "Subclusters",
    y = "Cells",
    fill = "Cell states"
  ) +
  theme_classic() 
dev.off()

## compute markers for the clusters
Idents(hb) <- hb$seurat_clusters
DefaultAssay(hb) <- "RNA"
hb_markers <- FindAllMarkers(hb, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)

# filter markers for adjusted p value < 0.05
hb_markers <- subset(hb_markers, p_val_adj < 0.05 & avg_log2FC > 1)
# add cell type name to markers 
hb_markers$cell_type_updated <- ifelse(hb_markers$cluster == 0, "Hindbrain 1 (late)", 
                                 ifelse(hb_markers$cluster == 1, "Hindbrain", "Hindbrain 2 (early)"))
write.csv(hb_markers, file="tables/scRNAseq/TOME/Hindbrain_subclustered_markers.csv")

## Do heatmap of top 10 markers per cluster
top30 <- hb_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

hb <- ScaleData(hb, verbose = FALSE)
pdf("images/scRNAseq/TOME/Hindbrain_subclustered_top30_markers_heatmap.pdf", width=10, height=14)
DoHeatmap(hb, features = top30$gene, group.by = "seurat_clusters", raster=T, size=3, angle=90)
dev.off()

fb_mb <- subset(integrated, subset = cell_type %in% c("Forebrain/midbrain"))
# subcluster processing
DefaultAssay(fb_mb)
fb_mb <- ScaleData(fb_mb, verbose = FALSE)
fb_mb <- RunPCA(fb_mb, npcs = 50, verbose = FALSE)
fb_mb <- FindNeighbors(fb_mb, dims = 1:50)
fb_mb <- FindClusters(fb_mb, resolution = 0.05)

plot_data <- fb_mb@meta.data
plot_data_summary <- plot_data %>%
  group_by(seurat_clusters, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(seurat_clusters) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  arrange(seurat_clusters, desc(percentage))

# plot
pdf("images/scRNAseq/TOME/stacked_Forebrain_midbrain_subclustered.pdf", width=8, height=6)
ggplot(plot_data_summary, aes(x = factor(seurat_clusters), y = count, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors) + 
  labs(
    x = "Subclusters",
    y = "Cells",
    fill = "Cell states"
  )+theme_classic()
dev.off()

## compute markers for the clusters
Idents(fb_mb) <- fb_mb$seurat_clusters
DefaultAssay(fb_mb) <- "RNA"
fb_mb_markers <- FindAllMarkers(fb_mb, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
# filter markers for adjusted p value < 0.05
fb_mb_markers <- subset(fb_mb_markers, p_val_adj < 0.05 & avg_log2FC > 1)
fb_mb_markers$cell_type_updated <- ifelse(fb_mb_markers$cluster == 0, "Midbrain", 
                                 ifelse(fb_mb_markers$cluster == 1, "Forebrain", "Unassigned"))
write.csv(fb_mb_markers, file="tables/scRNAseq/TOME/Forebrain_midbrain_subclustered_markers.csv")
## Do heatmap of top 10 markers per cluster
top30 <- fb_mb_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
DefaultAssay(fb_mb) <- "RNA"
fb_mb <- ScaleData(fb_mb, verbose = FALSE)

pdf("images/scRNAseq/TOME/Forebrain_midbrain_subclustered_top30_markers_heatmap.pdf", width=10, height=12)
DoHeatmap(fb_mb, features = top30$gene, group.by = "seurat_clusters", raster=T, size=3, angle=90)
dev.off()

DefaultAssay(integrated) <- "RNA"
integrated <- ScaleData(integrated, verbose = FALSE)

# from fb_mb cells rename the cell type to Forebrain/midbrain_0, 1,2
fb_mb$cell_type <- paste0("Forebrain/midbrain_", fb_mb$seurat_clusters)
fb_mb_vec <- setNames(fb_mb@meta.data$cell_type, Cells(fb_mb))
integrated <- AddMetaData(integrated, fb_mb_vec, col.name = "fb_mb_subcluster")

hb$cell_type <- paste0(hb$cell_type, "_", hb$seurat_clusters) 
hb_vec <- setNames(hb@meta.data$cell_type, Cells(hb))
integrated <- AddMetaData(integrated, hb_vec, col.name = "hb_subcluster")

integrated@meta.data$cell_type_updated <- ifelse(
  integrated@meta.data$cell_type == "Forebrain/midbrain",
  integrated@meta.data$fb_mb_subcluster,
  ifelse(
  integrated@meta.data$cell_type == "Hindbrain",
  integrated@meta.data$hb_subcluster,
integrated@meta.data$cell_type))
unique(integrated$cell_type_updated)

pdf('images/scRNAseq/TOME/UMAP_TOME_E725_E105_celltype_updated.pdf', width=12)
DimPlot(integrated, reduction = "umap", group.by='cell_type_updated', label=T,raster=F, cols=cell_type_updated_colors)
dev.off()

# rename annotations:
anno <- c(
           "Forebrain/midbrain_0" = "Midbrain",
           "Forebrain/midbrain_1" = "Forebrain",
           "Forebrain/midbrain_2" = "Unassigned")

# rename cell_type_updated to more descriptive names
integrated$cell_type_updated <- plyr::mapvalues(integrated$cell_type_updated, from = names(anno), to = anno)
# rename same cells in integrated 
integrated_vec <- setNames(integrated@meta.data$cell_type_updated, Cells(integrated))
integrated <- AddMetaData(integrated, integrated_vec, col.name = "cell_type_updated")
# take over cell_type in cell_type_updated 
integrated$cell_type_updated <- ifelse(is.na(integrated$cell_type_updated), integrated$cell_type, integrated$cell_type_updated)
integrated$cell_type_updated <- ifelse(integrated$cell_type == "Hindbrain", integrated$hb_subcluster, integrated$cell_type_updated)

# rename subclusters 
integrated$cell_type_updated <- ifelse(integrated$cell_type_updated == "Hindbrain_1", "Hindbrain", as.character(integrated$cell_type_updated))
integrated$cell_type_updated <- ifelse(integrated$cell_type_updated == "Unassigned 1", "Unassigned", as.character(integrated$cell_type_updated))
integrated$cell_type_updated <- ifelse(integrated$cell_type_updated == "Hindbrain_0", "Hindbrain 1", as.character(integrated$cell_type_updated)) # late Hindbrain
integrated$cell_type_updated <- ifelse(integrated$cell_type_updated == "Hindbrain_2", "Hindbrain 2", as.character(integrated$cell_type_updated)) # early Hindbrain

saveRDS(integrated, 'data/TOME/TOME_E7.25_E10.5_rpca.rds')

#########################################################################################################################################

## load gastruloids
Asmb <- readRDS('data/scRNAseq/Asmb.rds')
load(file="data/TLS/TLS_norm.Robj")
integrated <- readRDS('data/TOME/TOME_E7.25_E10.5_rpca.rds')
Asmb <- subset(Asmb, condition == "Hypoxic")

## Find anchors between reference and each query dataset (HAP, TLS_norm)
anchors_Asmb <- FindTransferAnchors(reference = integrated, query = Asmb, dims = 1:50, reference.reduction = "pca")  
anchors_TLS_norm <- FindTransferAnchors(reference = integrated, query = TLS_norm, dims = 1:50, reference.reduction = "pca")

Asmb_projected <- MapQuery(anchorset = anchors_Asmb, reference = integrated, query = Asmb, refdata = list(cell_type = "cell_type_updated", day = "day"), reference.reduction = "pca", reduction.model = "umap")
TLS_projected <- MapQuery(anchorset = anchors_TLS_norm, reference = integrated, query = TLS_norm, refdata = list(cell_type = "cell_type_updated", day = "day"), reference.reduction = "pca", reduction.model = "umap")

saveRDS(Asmb_projected, 'data/scRNAseq/Asmb_projected_Hypoxic_TOME_E7.25_E10.5_rpca.rds')
saveRDS(TLS_projected, 'data/scRNAseq/TLS_projected_Hypoxic_TOME_E7.25_E10.5_rpca.rds')

## Transfer day labels from reference to queries
Asmb_day <- Asmb_projected@meta.data[,c('condition','predicted.day','predicted.day.score','predicted.cell_type','predicted.cell_type.score')]
TLS_day <- TLS_projected@meta.data[,c('condition','predicted.day','predicted.day.score','predicted.cell_type','predicted.cell_type.score')]

day <- rbind(Asmb_day, TLS_day)
day$day.score <- ifelse(day$predicted.day.score > 0.6, '>0.6', '<=0.6')
day$predicted.day <- factor(day$predicted.day, levels=c("E7.25","E7.5","E7.75", "E8","E8.25","E8.5a", "E8.5b", "E9.5","E10.5"))

df <- rbind(Asmb_day, TLS_day)
df$day.score <- ifelse(df$predicted.day.score > 0.6, '>0.6', '<=0.6')
df$celltype.score <- ifelse(df$predicted.cell_type.score > 0.6, '>0.6', '<=0.6')
df$predicted.day <- factor(df$predicted.day, levels=c("E7.25","E7.5","E7.75", "E8","E8.25","E8.5a", "E8.5b", "E9.5","E10.5"))

## plot day assignment
pdf('images/scRNAseq/TOME/Day_prediction_TOME_E7.25_E10.5_rpca_days.pdf', width=18, height=5)
ggplot(day, aes(x=predicted.day, fill=day.score)) + geom_bar() + theme_classic() + facet_grid(~condition)
dev.off()


## plot cell type assignment
pdf('images/scRNAseq/TOME/CellType_prediction_TOME_E7.25_E10.5_rpca.pdf', width=10)
ggplot(df, aes(x=predicted.cell_type, fill=celltype.score)) + geom_bar() + theme_classic() + facet_grid(condition~.) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('images/scRNAseq/TOME/CellType_prediction_per_day_TOME_E7.25_E10.5_rpca.pdf', width=28, height=10)
ggplot(df, aes(x=predicted.cell_type, fill=celltype.score)) + geom_bar() + theme_classic() + facet_grid(predicted.day~condition) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# filter
combined <- merge(Asmb_projected, TLS_projected)
df <- combined@meta.data

# Cell type prediction score
score_cutoff <- quantile(df$predicted.cell_type.score, 0.2, na.rm = TRUE)
pdf('images/scRNAseq/TOME/Distr_PredictionScore_celltype_20.pdf', width=8, height=5)
ggplot(df, aes(x=predicted.cell_type.score)) +
  geom_histogram(bins=50, fill="#1f77b4", color="black") +
  geom_vline(xintercept = score_cutoff, color="red", linetype="dashed", size=1) +
  theme_classic() +
  labs(x="Predicted Cell Type Score", y="Cell Count", title="Distribution of Cell Type Prediction Scores") +
  annotate("text", x=score_cutoff, y=Inf, label=paste0("20th percentile: ", round(score_cutoff, 3)), vjust=2, hjust=0, color="red")
dev.off()

# Day prediction score
score_cutoff <- quantile(df$predicted.day.score, 0.2, na.rm = TRUE)
pdf('images/scRNAseq/TOME/Distr_PredictionScore_day_20.pdf', width=8, height=5)
ggplot(df, aes(x=predicted.day.score)) +
  geom_histogram(bins=50, fill="#2ca02c", color="black") +
  geom_vline(xintercept = score_cutoff, color="red", linetype="dashed", size=1) +
  theme_classic() +
  labs(x="Predicted Day Score", y="Cell Count", title="Distribution of Day Prediction Scores") +
  annotate("text", x=score_cutoff, y=Inf, label=paste0("20th percentile: ", round(score_cutoff, 3)), vjust=2, hjust=0, color="red")
dev.off()

## apply prediction score filter 
## filter out all the cells with the predicted.celltype_update.score < 20th percentile
percentile_celltype <- quantile(df$predicted.cell_type.score, 0.2, na.rm = TRUE)
percentile_day <- quantile(df$predicted.day.score, 0.2, na.rm = TRUE)

Asmb_projected_f <- subset(Asmb_projected, subset = predicted.cell_type.score >= percentile_celltype & predicted.day.score >= percentile_day)
TLS_projected_f <- subset(TLS_projected, subset = predicted.cell_type.score >= percentile_celltype  & predicted.day.score >= percentile_day)

combined <- merge(Asmb_projected_f, TLS_projected_f)
df <- combined@meta.data

orig_df <- df  
ct10_by_cond <- orig_df %>%
  group_by(condition, predicted.cell_type) %>%
  tally(name = "n") %>%
  filter(n >= 10) %>%
  arrange(condition, desc(n))

print(ct10_by_cond, n=41)

allowed_pairs <- ct10_by_cond %>%
  mutate(key = paste(condition, predicted.cell_type, sep = "||")) %>%
  pull(key)

meta_asmb <- Asmb_projected_f@meta.data
meta_asmb$key <- paste(meta_asmb$condition, meta_asmb$predicted.cell_type, sep = "||")
keep_asmb <- rownames(meta_asmb)[meta_asmb$key %in% allowed_pairs]
Asmb_projected_f <- subset(Asmb_projected_f, cells = keep_asmb)

# subset TLS_projected_f similarly
meta_tls <- TLS_projected_f@meta.data
meta_tls$key <- paste(meta_tls$condition, meta_tls$predicted.cell_type, sep = "||")
keep_tls <- rownames(meta_tls)[meta_tls$key %in% allowed_pairs]
TLS_projected_f <- subset(TLS_projected_f, cells = keep_tls)

table(Asmb_projected_f$condition, Asmb_projected_f$predicted.cell_type)
table(TLS_projected_f$condition, TLS_projected_f$predicted.cell_type)

filtered <- merge(Asmb_projected_f, TLS_projected_f)
filtered_meta <- filtered@meta.data

# DotPlot to show cell fraction per predicted stage 
dotplot_data <- filtered_meta %>%
  dplyr::group_by(condition, predicted.day) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%  # normalization within each condition
  dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) 
print(dotplot_data)

## condition pairs for comparisons
condition_pairs <- list(
  c("TLS", "Hypoxic")
)

all_stages <- c("E7.25","E7.5","E7.75", "E8","E8.25","E8.5a", "E8.5b", "E9.5","E10.5")

# dot plots for each condition pair
for (pair in condition_pairs) {
  
  dotplot_subset <- filtered_meta %>%
    filter(condition %in% pair) %>%
    dplyr::group_by(condition, predicted.day) %>%
    dplyr::summarise(cell_count = n(), .groups = "drop") %>%
    dplyr::group_by(condition) %>%  
    dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) %>%
    ungroup()  
  dotplot_subset <- dotplot_subset %>%
    tidyr::complete(condition = pair, predicted.day = all_stages, fill = list(cell_count = NA, percentage = NA))
  dotplot_subset$predicted.day <- factor(dotplot_subset$predicted.day, levels = all_stages, ordered = TRUE)
  
  plot <- ggplot(dotplot_subset, aes(x = predicted.day, y = condition, size = percentage)) +
    geom_point(alpha = 0.8, na.rm = TRUE) + 
    scale_size_continuous(range = c(2, 10)) +  # No dots for 0% values
    scale_x_discrete(expand = c(0.1, 0.1)) +
     coord_cartesian(clip = "off") +
    labs(x = "Predicted stage (E)", y = "", size = "Cell fraction (%)", title = "") +
    theme_minimal() +
    theme_classic() +
    theme(
      # rotate x-axis labels
      axis.text.x = element_text(size = 14, color = "black", angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),  
      axis.title = element_text(size = 16),  
      legend.text = element_text(size = 12),  
      legend.title = element_text(size = 14), aspect.ratio = 3/5
    )
  filename <- paste0("images/scRNAseq/TOME/TOME_20_dotplot_stage_", pair[1], "_vs_", pair[2], ".pdf")
  ggsave(filename, width = 7, height = 3, plot)
  print(pair)
}

filtered_meta$predicted.days <- ifelse(filtered_meta$predicted.day %in% c("E8.5a", "E8.5b"), "E8.5", as.character(filtered_meta$predicted.day))

# DotPlot to show cell fraction per predicted stage 
dotplot_data <- filtered_meta %>%
  dplyr::group_by(condition, predicted.days) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%  # normalization within each condition
  dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) 
print(dotplot_data)

all_stages <- c("E7.25","E7.5","E7.75", "E8","E8.25","E8.5", "E9.5","E10.5")

# dot plots for each condition pair
for (pair in condition_pairs) {

  dotplot_subset <- filtered_meta %>%
    filter(condition %in% pair) %>%
    dplyr::group_by(condition, predicted.days) %>%
    dplyr::summarise(cell_count = n(), .groups = "drop") %>%
    dplyr::group_by(condition) %>%  
    dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) %>%
    ungroup()  
  dotplot_subset <- dotplot_subset %>%
    tidyr::complete(condition = pair, predicted.days = all_stages, fill = list(cell_count = NA, percentage = NA))
  dotplot_subset$predicted.days <- factor(dotplot_subset$predicted.days, levels = all_stages, ordered = TRUE)
  
  plot <- ggplot(dotplot_subset, aes(x = predicted.days, y = condition, size = percentage)) +
    geom_point(alpha = 0.8, na.rm = TRUE) + 
    scale_size_continuous(range = c(2, 10)) +  # No dots for 0% values
    scale_x_discrete(expand = c(0.1, 0.1)) +
     coord_cartesian(clip = "off") +
    labs(x = "Predicted stage (E)", y = "", size = "Cell fraction (%)", title = "") +
    theme_minimal() +
    theme_classic() +
    theme(
      # rotate x-axis labels for better readability
      axis.text.x = element_text(size = 14, color = "black", angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),  
      axis.title = element_text(size = 16),  
      legend.text = element_text(size = 12),  
      legend.title = element_text(size = 14), aspect.ratio = 3/5
    )
  filename <- paste0("images/scRNAseq/TOME/TOME_20_dotplot_stage_8.5fused_", pair[1], "_vs_", pair[2], ".pdf")
  ggsave(filename, width = 7, height = 3, plot)
  print(pair)
}

#############################################################################################################################
# genrate UMAP of TOME using only cell states that are present in filtered

integrated <- readRDS('data/TOME/TOME_E7.25_E10.5_rpca.rds')

a <- unique(Asmb_projected_f$predicted.cell_type)
b <- unique(TLS_projected_f$predicted.cell_type)
ct10 <- union(a, b)
ct10 <- as.character(ct10)

integrated_f <- subset(integrated, subset = cell_type_updated %in% ct10)
unique(integrated_f$cell_type_updated)
DefaultAssay(integrated_f)

# reprocess filtered_OMG
integrated_f <- ScaleData(integrated_f, verbose = FALSE)
integrated_f <- RunPCA(integrated_f, npcs = 50, verbose = FALSE)
integrated_f <- RunUMAP(integrated_f, dims = 1:50, return.model = T)

source("scripts/TOME_colors.r")

unique(integrated_f$cell_type_updated)
pdf('images/scRNAseq/TOME/UMAP_filtered_TOME_E725_E105_celltype.pdf', width=10)
DimPlot(integrated_f, reduction = "umap", group.by='cell_type_updated', raster=F, cols=cell_type_updated_colors)
dev.off()

unique(integrated_f$day)
pdf('images/scRNAseq/TOME/UMAP_filtered_TOME_E725_E105_day.pdf', width=10)
DimPlot(integrated_f, reduction = "umap", group.by='day', raster=F, cols=day_colors, shuffle=TRUE)
dev.off()


# project filtered TLS and HAP onto UMAP of integrated_f

Asmb_projected_f <- NormalizeData(Asmb_projected_f, verbose = FALSE)
Asmb_projected_f <- FindVariableFeatures(Asmb_projected_f, verbose = FALSE)
Asmb_projected_f <- ScaleData(Asmb_projected_f, verbose = FALSE)
Asmb_projected_f <- RunPCA(Asmb_projected_f, npcs = 50, verbose = FALSE)

TLS_projected_f <- NormalizeData(TLS_projected_f, verbose = FALSE)  
TLS_projected_f <- FindVariableFeatures(TLS_projected_f, verbose = FALSE)
TLS_projected_f <- ScaleData(TLS_projected_f, verbose = FALSE)
TLS_projected_f <- RunPCA(TLS_projected_f, npcs = 50, verbose = FALSE)

anchors_Asmb <- FindTransferAnchors(reference = integrated_f, query = Asmb_projected_f, dims = 1:50, reference.reduction = "pca")  
anchors_TLS_normm <- FindTransferAnchors(reference = integrated_f, query = TLS_projected_f, dims = 1:50, reference.reduction = "pca")
Asmb_projected_f <- MapQuery(anchorset = anchors_Asmb, reference = integrated_f, query = Asmb_projected_f,  reference.reduction = "pca", reduction.model = "umap")
TLS_projected_f <- MapQuery(anchorset = anchors_TLS_normm, reference = integrated_f, query = TLS_projected_f, reference.reduction = "pca", reduction.model = "umap")

# check of Asmb projected umap 
umap_data_new <- as.data.frame(Asmb_projected_f[["ref.umap"]]@cell.embeddings)
head(umap_data_new)
umap_data_WT <- as.data.frame(integrated_f[["umap"]]@cell.embeddings)
head(umap_data_WT) ##  umap_1 
umap_data_old <- as.data.frame(TLS_projected_f[["ref.umap"]]@cell.embeddings)
head(umap_data_old)

umap_data_old$condition <- TLS_projected_f@meta.data[row.names(umap_data_old),"condition"]
umap_data_new$condition <- Asmb_projected_f@meta.data[row.names(umap_data_new),"condition"]
umap_data_old$cell_state_annot <- TLS_projected_f@meta.data[row.names(umap_data_old),"predicted.cell_type"]
umap_data_new$cell_state_annot <- Asmb_projected_f@meta.data[row.names(umap_data_new),"predicted.cell_type"]
umap_data_old$day <- TLS_projected_f@meta.data[row.names(umap_data_old),"predicted.day"]
umap_data_new$day <- Asmb_projected_f@meta.data[row.names(umap_data_new),"predicted.day"]

# rename refUMAP columns to match integrated UMAP columns
colnames(umap_data_new)[1:2] <- c("refumap_1", "refumap_2")

umap_combined <- bind_rows(
  umap_data_old %>% mutate(source = "TLS"),
  umap_data_new %>% mutate(source = "Asmb")
)
umap_combined$condition <- factor(umap_combined$condition)

umap_data_WT_rep <- umap_data_WT[rep(1:nrow(umap_data_WT), times = length(unique(umap_combined$condition))), ]
umap_data_WT_rep$condition <- rep(unique(umap_combined$condition), each = nrow(umap_data_WT))

label_positions <- umap_combined %>%
  group_by(cell_state_annot) %>%
  summarise(refumap_1 = mean(refumap_1), refumap_2 = mean(refumap_2)) %>%
  mutate(cell_state_id = str_extract(cell_state_annot, "\\d+"))  # Extract only numbers

# Define unique conditions
conditions <- unique(umap_combined$condition)

umap_combined$cell_state_annot <- factor(
  umap_combined$cell_state_annot,
  levels = names(cell_type_updated_colors)
)
umap_combined$cell_type_num <- unname(
  cell_type_ids[ as.character(umap_combined$cell_state_annot) ]
)
umap_combined$cell_type_num <- factor(
  umap_combined$cell_type_num,
  levels = unname(cell_type_ids)   
)
umap_combined$num <- str_extract(umap_combined$cell_type_num, "\\d+")


# Loop through each condition and save a separate plot
for (cond in conditions) {
  
  # Filter data for the current condition
  umap_subset <- subset(umap_combined, condition == cond)
  
    centroids <- umap_subset %>%
    dplyr::group_by(cell_type_num, num) %>%
    dplyr::summarise(refumap_1 = median(refumap_1), refumap_2 = median(refumap_2), .groups = "drop")

  # Generate plot
  single_plot <- ggplot() +
    rasterise(
      geom_point(data = umap_data_WT_rep, aes(x = umap_1, y = umap_2), 
                 size = 3, colour = "grey82", alpha = 0.5),
      dpi = 300
    ) +
    geom_point(data = umap_subset, aes(x = refumap_1, y = refumap_2, color = cell_type_num), 
               size = 3) +
    scale_color_manual(name = "", values = cell_type_final_colors) +
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
    geom_text(data = centroids, aes(x = refumap_1, y = refumap_2, label = num),
              size = 6, fontface = "bold", color = "black") +
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
  
  # Save plot
  ggsave(paste0("images/scRNAseq/TOME/UMAP_", cond, "_projected_f20_TOME.pdf"),
         single_plot, width = 10, height = 7)

  print(cond)
  print(dim(umap_subset))
}


head(umap_combined)
umap_combined$day_f <- ifelse(umap_combined$day %in% c("E8.5a", "E8.5b"), "E8.5", as.character(umap_combined$day))
umap_combined$day_f <- factor(umap_combined$day_f, levels=c("E7.25","E7.5","E7.75", "E8","E8.25","E8.5", "E9.5","E10.5"))
source("scripts/TOME_colors.r")

# Loop through each condition and save a separate plot
for (cond in conditions) {
  
  # Filter data for the current condition
  umap_subset <- subset(umap_combined, condition == cond)
  
  # Generate plot
  single_plot <- ggplot() +
    rasterise(
      geom_point(data = umap_data_WT_rep, aes(x = umap_1, y = umap_2), 
                 size = 3, colour = "grey82", alpha = 0.5),
      dpi = 300
    ) +
    geom_point(data = umap_subset, aes(x = refumap_1, y = refumap_2, color = day_f), 
               size = 3) +
    scale_color_manual(name = "", values = day_colors_f) +
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
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
  
  # Save plot
  ggsave(paste0("images/scRNAseq/TOME/UMAP_", cond, "_projected_day_f20_TOME.pdf"),
         single_plot, width = 10, height = 7)

  print(cond)
  print(dim(umap_subset))
}

unique(integrated_f$cell_type_updated)

umap_data_new$predicted.cell_type.score <- Asmb_projected_f@meta.data[row.names(umap_data_new),"predicted.cell_type.score"]
umap_data_old$predicted.cell_type.score <- TLS_projected_f@meta.data[row.names(umap_data_old),"predicted.cell_type.score"]

umap_combined <- bind_rows(
  umap_data_old %>% mutate(source = "TLS"),
  umap_data_new %>% mutate(source = "Asmb")
)
umap_combined$condition <- factor(umap_combined$condition)



# plot umap colored by cell_type.score, lowest on top 
for(cond in conditions) {
  
  # Filter data for the current condition
  umap_subset <- subset(umap_combined, condition == cond)
  
  # Generate plot
  single_plot <- ggplot() +
  rasterise(
      geom_point(data = umap_data_WT_rep, aes(x = umap_1, y = umap_2), 
                 size = 3, colour = "grey82", alpha = 0.5),
      dpi = 300
    ) +
    geom_point(data = umap_subset %>% arrange(predicted.cell_type.score), aes(x = refumap_1, y = refumap_2, color = predicted.cell_type.score), 
               size = 3) +
    scale_color_viridis_c(name = "Cell type\nprediction score", option = "magma", direction = -1) +
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
    guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    coord_fixed(ratio = 1) 
  
  # Save plot
  ggsave(paste0("images/scRNAseq/TOME/UMAP_", cond, "_projected_f20_TOME_celltype_score.pdf"),
         single_plot, width = 10, height = 7)
  
  print(cond)
  print(dim(umap_subset))
} 

umap_data_WT$cell_type <-  integrated_f@meta.data[row.names(umap_data_WT),"cell_type_updated"]
umap_data_WT$cell_type_num <- unname(
  cell_type_ids[ as.character(umap_data_WT$cell_type) ]
)
umap_data_WT$cell_type_num <- factor(
  umap_data_WT$cell_type_num,
  levels = unname(cell_type_ids)   # "Neuromesodermal progenitors (1)", ..., "Unassigned 2 (20)"
)
umap_data_WT$num <- str_extract(umap_data_WT$cell_type_num, "\\d+")
centroids <- umap_data_WT %>%
dplyr::group_by(cell_type_num, num) %>%
dplyr::summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")


p<- ggplot() +
    geom_point(data = umap_data_WT, aes(x = umap_1, y = umap_2, color = cell_type_num), 
               size = 3) +
    scale_color_manual(name = "", values = cell_type_final_colors) +
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
  
  # Save plot
ggsave(paste0("images/scRNAseq/TOME/UMAP_ref_nums_f20_TOME.pdf"), p, width = 10, height = 7)

# [1] "TLS"
# [1] 2438    9
# [1] "Hypoxic"
# [1] 1379    9


umap_data_WT$day <- integrated_f@meta.data[row.names(umap_data_WT),"day"]
umap_data_WT$day_f <- ifelse(umap_data_WT$day %in% c("E8.5a", "E8.5b"), "E8.5", as.character(umap_data_WT$day))
umap_data_WT$day_f <- factor(umap_data_WT$day_f, levels=c("E7.25","E7.5","E7.75", "E8","E8.25","E8.5", "E9.5","E10.5"))

umap_data_WT_shuf <- umap_data_WT[sample(nrow(umap_data_WT)), ]

p<- ggplot() +
    geom_point(data = umap_data_WT_shuf, aes(x = umap_1, y = umap_2, color = day_f), 
               size = 3) +
    scale_color_manual(name = "", values = day_colors_f) +
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
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
  
  # Save plot
ggsave(paste0("images/scRNAseq/TOME/UMAP_ref_day_nums_f20_TOME.pdf"), p, width = 10, height = 7)


# stacked bar plots
calculate_cell_state_percentage <- function(seurat_object) {
  seurat_object@meta.data %>%
    dplyr::group_by(condition, predicted.cell_type) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(percentage = count / sum(count) * 100) %>%
    dplyr::arrange(condition, dplyr::desc(percentage))
}

# Compute cell state percentages for both objects
tls_plot_data <- calculate_cell_state_percentage(TLS_projected_f)
asmb_plot_data <- calculate_cell_state_percentage(Asmb_projected_f)
plot_data <- bind_rows(tls_plot_data, asmb_plot_data)

# Ensure ordered levels for plotting
# need them to be ordered as in cell_type_colors
all_levels     <- names(cell_type_updated_colors)
present_levels <- intersect(all_levels, unique(plot_data$predicted.cell_type))
plot_data$predicted.cell_type <- factor(plot_data$predicted.cell_type,
                                              levels = present_levels)

# rename plot_data Hypoxic to HAP
plot_data$condition <- ifelse(plot_data$condition == "Hypoxic", "HAP", as.character(plot_data$condition))
plot_data$condition <- factor(plot_data$condition, levels = c("TLS", "HAP"))
# Plot
pdf("images/scRNAseq/TOME/stacked_TOME_f20_predicted_states_conditions.pdf", width=10)
ggplot(plot_data, aes(x = factor(condition), y = percentage, fill = predicted.cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = cell_type_updated_colors,
    breaks = present_levels   # ensures legend order matches
  ) +
  labs(
    x = "",
    y = "Percentage of cells",
    fill = "Cell states"
  ) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  )

dev.off()

combined <- merge(Asmb_projected_f, TLS_projected_f)
dim(combined)
df <- combined@meta.data
df$predicted.days <- ifelse(df$predicted.day %in% c("E8.5a", "E8.5b"), "E8.5", as.character(df$predicted.day))

# DotPlot to show cell fraction per predicted stage 
dotplot_data <- df %>%
  dplyr::group_by(condition, predicted.days) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%  # normalization within each condition
  dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) 
print(dotplot_data)

all_stages <- c("E7.25","E7.5","E7.75", "E8","E8.25","E8.5", "E9.5","E10.5")
# rename Hypoxic condition to HAP
df$condition <- ifelse(df$condition == "Hypoxic", "HAP", as.character(df$condition))
df$condition <- factor(df$condition, levels = c("TLS", "HAP"))
condition_pairs <- list(
  c("TLS", "HAP")
)
# dot plots for each condition pair
for (pair in condition_pairs) {
  
  dotplot_subset <- df %>%
    filter(condition %in% pair) %>%
    dplyr::group_by(condition, predicted.days) %>%
    dplyr::summarise(cell_count = n(), .groups = "drop") %>%
    dplyr::group_by(condition) %>%  
    dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) %>%
    ungroup()  
  dotplot_subset <- dotplot_subset %>%
    tidyr::complete(condition = pair, predicted.days = all_stages, fill = list(cell_count = NA, percentage = NA))
  dotplot_subset$predicted.days <- factor(dotplot_subset$predicted.days, levels = all_stages, ordered = TRUE)
  
  plot <- ggplot(dotplot_subset, aes(x = predicted.days, y = condition, size = percentage)) +
    geom_point(alpha = 0.8, na.rm = TRUE) + 
    scale_size_continuous(range = c(2, 10)) +  # No dots for 0% values
    scale_x_discrete(expand = c(0.1, 0.1)) +
     coord_cartesian(clip = "off") +
    labs(x = "Predicted stage (E)", y = "", size = "Cell fraction (%)", title = "") +
    theme_minimal() +
    theme_classic() +
    theme(
      # rotate x-axis labels for better readability
      axis.text.x = element_text(size = 14, color = "black", angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),  
      axis.title = element_text(size = 16),  
      legend.text = element_text(size = 12),  
      legend.title = element_text(size = 14), aspect.ratio = 3/5
    )
  filename <- paste0("images/scRNAseq/TOME/TOME_20_dotplot_stage_8.5fused_", pair[1], "_vs_", pair[2], ".pdf")
  ggsave(filename, width = 7, height = 3, plot)
  print(pair)
}


# Staging per cell type 
day_scores <- GetAssayData(Asmb_projected_f, assay = "prediction.score.day", layer = "data")
haps <- subset(Asmb_projected_f, condition == "Hypoxic")
hapa_cells <- colnames(haps)
day_scores <- day_scores[, hapa_cells]

w <- parse_number(rownames(day_scores)) 
# Weighted sum (expected somite count) per cell
num   <- colSums(sweep(day_scores, 1, w, `*`))      # sum (score * somite)                
expected_day <- num 

df_w <- data.frame(
  cell = colnames(day_scores),
  expected_day = expected_day,
  row.names = colnames(day_scores)
)

meta <- haps@meta.data[colnames(day_scores), , drop = FALSE]  
df_w$cell_type <- meta$predicted.cell_type
df_w$cell_type <- factor(df_w$cell_type,
                         levels = names(cell_type_ids))

yticks <- c(8.0, 8.5, 9.0, 9.5, 10.0)

pdf("images/scRNAseq/TOME/TOME_20_boxplot_day_weighted_sum_HAPs_per_celltype.pdf",
    width = 12, height = 8)
ggplot(df_w, aes(x = cell_type, y = expected_day,  fill = cell_type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(y = "Predicted stage", x = "") +
   scale_fill_manual(values = cell_type_updated_colors) +
   scale_y_continuous(
    breaks = yticks,
    limits = c(8, 9.5),
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 60), 
  )
dev.off()


# plot markers for filtered TOME reference 
Idents(integrated_f) <- "cell_type_updated"
DefaultAssay(integrated_f) <- "RNA"
markers <- FindAllMarkers(integrated_f, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
# filter markers with adjusted p value < 0.05
markers <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 1, ]
write.csv(markers, file="tables/scRNAseq/TOME/TOME_filtered_markers.csv")
write.xlsx(markers, file="tables/scRNAseq/TOME/TOME_filtered_markers.xlsx")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# z-scores 
integrated_f <- ScaleData(
  integrated_f,
  assay   = "RNA",
  features = top10$gene,
  verbose = FALSE
)
top10 <- top10 %>%
  mutate(cluster = factor(cluster, levels = names(cell_type_ids))) %>%
  arrange(cluster, desc(avg_log2FC))
# keep cluster blocks separated 
feature_list <- split(top10$gene, top10$cluster)

Idents(integrated_f) <- factor(
  Idents(integrated_f),
  levels = names(cell_type_ids)  #
)

#heatmap
pdf("images/scRNAseq/TOME/Heatmap_marker_expression_TOME_filtered_top10.pdf", width=20, height=18)  
DoHeatmap(integrated_f, features = top10$gene, size=4, group.colors=cell_type_updated_colors, raster=TRUE) + 
  theme(axis.text.y = element_text(size = 5))
dev.off()

## save objects
saveRDS(integrated_f, file="data/scRNAseq/TOME_filtered_20_TOME_E725_E105.rds")
saveRDS(Asmb_projected_f, file="data/scRNAseq/Asmb_projected_filtered_20_TOME_E725_E105.rds")
saveRDS(TLS_projected_f, file="data/scRNAseq/TLS_projected_filtered_20_TOME_E725_E105.rds")
saveRDS(integrated, file="data/TOME/TOME_E7.25_E10.5_rpca.rds")

TLS_projected_f <- readRDS("data/scRNAseq/TLS_projected_filtered_20_TOME_E725_E105.rds")
Asmb_projected_f <- readRDS("data/scRNAseq/Asmb_projected_filtered_20_TOME_E725_E105.rds")
integrated_f <- readRDS("data/scRNAseq/TOME_filtered_20_TOME_E725_E105.rds")