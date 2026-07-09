# Gastruloid scRNA-seq processing and OMG-based cell-state, stage and somite annotation
#
# Input : Parse filtered matrices for four conditions
#         (Normoxic / Hypoxic / Hypo+XAV / HIF1A-KO); the annotated OMG
#         reference (data/OMG/OMG_E8_E9.75_merge.rds, from 02) and the
#         normoxic TLS reference.
#
# Steps : 1. Load per-condition matrices, merge into one Seurat object
#         2. QC filtering (nFeature / nCount / percent.mt) + scDblFinder;
#            writes data/scRNAseq/Asmb.rds 
#         3. Normalize, scale, PCA, UMAP
#         4. Label transfer from OMG (cell state, day, somite count) via
#            FindTransferAnchors / MapQuery, filter by prediction score
#            (<20th percentile) and drop states with <10 cells
#         5. Project each condition onto the OMG reference UMAP
#         6. Cell-state composition, predicted stage /
#            somite-count summaries, marker dot plots and heatmaps
#            (neural, gut, somite lineages)
#
# Output: Seurat objects in data/scRNAseq/ (Asmb, Asmb_OMG, OMG_filtered,
#         TLS_OMG, Asmb_neural, ...), figures in images/scRNAseq/,
#         marker tables in tables/scRNAseq/OMG/

###############################################################################
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
library(scDblFinder)
library(SingleCellExperiment)
library(tidyr)
library(openxlsx)
set.seed(42)
source("scripts/OMG_colors.r")


###############################################################################
# Input matrices (Parse Evercode, per condition).
# Processed matrices are on GEO (GSE337584) as flat, condition-prefixed files:
#     NAP_count_matrix.mtx.gz  NAP_all_genes.csv.gz  NAP_cell_metadata.csv.gz
#     HAP_*  HAPX_*  HIF1AKO_*
# Condition-name mapping (GEO -> internal):
#     NAP = Normoxic | HAP = Hypoxic | HAPX = Hypo_XAV | HIF1AKO = HIF1AKO
# To run this script on the GEO files, place each condition's three files in
#   data/scRNAseq/filtered_matrices/output_combined/<condition>/DGE_filtered/
# renamed to count_matrix.mtx.gz / cell_metadata.csv.gz / all_genes.csv.gz
# (or edit `base` and `read_parse_sample()` to read the flat GEO names directly).
###############################################################################

base <- "data/scRNAseq/filtered_matrices/output_combined"

samples <- c(
  Normoxic = "Normoxic_aggregates___Tasos___ABK_Lab",
  Hypoxic  = "Hypoxic_aggregates___Tasos___ABK_Lab",
  Hypo_XAV = "Hypo_XAV_aggregates___Tasos___ABK_Lab",
  HIF1AKO  = "HIF1AKO_aggregates___Tasos___ABK_Lab"
)
read_parse_sample <- function(sample_dir, sample_name) {
  d <- file.path(sample_dir, "DGE_filtered")
  
  pick <- function(stem) {
    gz <- file.path(d, paste0(stem, ".gz"))
    if (file.exists(gz)) gz else file.path(d, stem)
  }
  
  mat <- ReadMtx(
    mtx            = pick("count_matrix.mtx"),
    cells          = pick("cell_metadata.csv"),
    features       = pick("all_genes.csv"),
    cell.column    = 1,
    feature.column = 2,
    cell.sep       = ",",      # <-- added
    feature.sep    = ",",      # <-- added
    skip.cell      = 1,
    skip.feature   = 1,
    mtx.transpose  = TRUE
  )
  
  meta <- read.csv(pick("cell_metadata.csv"), row.names = 1)
  
  so <- CreateSeuratObject(
    counts    = mat,
    meta.data = meta,
    project   = sample_name,
    min.cells    = 3,
    min.features = 200
  )
  so$condition <- sample_name
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
  so
}

so_list <- lapply(names(samples), function(nm) {
  message("Loading ", nm)
  read_parse_sample(file.path(base, samples[[nm]]), nm)
})
names(so_list) <- names(samples)

so_list
sapply(so_list, ncol)
sapply(so_list, function(x) median(x$percent.mt))


Asmb <- merge(
  so_list[[1]],
  y = so_list[-1],
  add.cell.ids = names(so_list),
  project = "Tasosoids"
)

Asmb
table(Asmb$condition)
table(Asmb$orig.ident)

Asmb <- JoinLayers(Asmb)
Asmb <- NormalizeData(Asmb)

Idents(Asmb) <- "condition"
# Visualize QC metrics as a violin plot
pdf("images/scRNAseq/Asmb_QC_violin.pdf", width = 8, height = 10)
Seurat::VlnPlot(Asmb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

# QC filtering
min.genes <- 300
max.genes <- 9000
min.UMI <- 700 
max.UMI <- 30000
max.mt    <- 5

Asmb <- subset(Asmb, subset = nFeature_RNA > min.genes & nFeature_RNA < max.genes & 
                          nCount_RNA < max.UMI & nCount_RNA > min.UMI & percent.mt < max.mt)

counts <- GetAssayData(Asmb, assay = "RNA", layer = "counts")  # v5-safe; avoids as.SCE quirks
set.seed(42)
sce <- scDblFinder(counts, samples = Asmb$condition)

# column order is preserved -> assign back by position
Asmb$scDblFinder.class <- sce$scDblFinder.class   # factor: "singlet" / "doublet"
Asmb$scDblFinder.score <- sce$scDblFinder.score

table(Asmb$condition, Asmb$scDblFinder.class)
Asmb <- subset(Asmb, subset = scDblFinder.class == "singlet")

pdf("images/scRNAseq/Asmb_QC_violin_subset.pdf", width = 8, height = 10)
Seurat::VlnPlot(Asmb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

Asmb <- FindVariableFeatures(Asmb, nfeatures = 2000)
all.genes <- rownames(Asmb)
Asmb <- ScaleData(Asmb, features = all.genes)
Asmb <- RunPCA(Asmb, npcs = 50, features = VariableFeatures(Asmb))

pdf("images/scRNAseq/Asmb_QC_elbow.pdf", width = 8, height = 4)
ElbowPlot(Asmb, ndims = 60)
dev.off()
Asmb <- RunUMAP(Asmb, dims = 1:30, return.model=TRUE)

pdf("images/scRNAseq/Asmb_QC_umap.pdf", width = 8, height = 4)
DimPlot(Asmb, group.by = "condition")
dev.off()

# featureplot for counts and mt.percent
pdf("images/scRNAseq/Asmb_QC_featureplot.pdf", width = 12, height = 4)
FeaturePlot(Asmb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# save Asmb object before mapping onto OMG
saveRDS(Asmb, file = "data/scRNAseq/Asmb.rds")
table(Asmb$condition)
OMG <- readRDS('data/OMG/OMG_E8_E9.75_merge.rds')
load(file="data/TLS/TLS_norm.Robj")

anchors_Asmb <- FindTransferAnchors(reference = OMG, query = Asmb, dims = 1:50, reference.reduction = "pca")  
anchors_TLS_norm <- FindTransferAnchors(reference = OMG, query = TLS_norm, dims = 1:50, reference.reduction = "pca")

Asmb <- MapQuery(anchorset = anchors_Asmb, reference = OMG, query = Asmb, refdata = list(day='day',somite_count='somite_count',major_trajectory='major_trajectory',celltype_update='celltype_update',celltype_updated='celltype_updated',neurons_sub_clustering='neurons_sub_clustering',neurons_sub_clustering='neurons_sub_clustering',lateral_plate_mesoderm_sub_clustering='lateral_plate_mesoderm_sub_clustering'), reference.reduction = "pca", reduction.model = "umap")
TLS_projected <- MapQuery(anchorset = anchors_TLS_norm, reference = OMG, query = TLS_norm, refdata = list(day='day',somite_count='somite_count',major_trajectory='major_trajectory',celltype_update='celltype_update',celltype_updated='celltype_updated',neurons_sub_clustering='neurons_sub_clustering',neurons_sub_clustering='neurons_sub_clustering',lateral_plate_mesoderm_sub_clustering='lateral_plate_mesoderm_sub_clustering'), reference.reduction = "pca", reduction.model = "umap")

combined <- merge(Asmb, TLS_projected)
df <- combined@meta.data

# Cell type prediction score
score_cutoff <- quantile(df$predicted.celltype_updated.score, 0.2, na.rm = TRUE)
pdf('images/scRNAseq/OMG/Distr_PredictionScore_celltype_updated_20.pdf', width=8, height=5)
ggplot(df, aes(x=predicted.celltype_updated.score)) +
  geom_histogram(bins=50, fill="#1f77b4", color="black") +
  geom_vline(xintercept = score_cutoff, color="red", linetype="dashed", size=1) +
  theme_classic() +
  labs(x="Predicted Cell Type Score", y="Cell Count", title="Distribution of Cell Type Prediction Scores") +
  annotate("text", x=score_cutoff, y=Inf, label=paste0("20th percentile: ", round(score_cutoff, 3)), vjust=2, hjust=0, color="red")
dev.off()

# Somite count prediction score
score_cutoff <- quantile(df$predicted.somite_count.score, 0.2, na.rm = TRUE)
pdf('images/scRNAseq/OMG/Distr_PredictionScore_somite_count_20.pdf', width=8, height=5)
ggplot(df, aes(x=predicted.somite_count.score)) +
  geom_histogram(bins=50, fill="#ff7f0e", color="black") +
  geom_vline(xintercept = score_cutoff, color="red", linetype="dashed", size=1) +
  theme_classic() +
  labs(x="Predicted Somite Count Score", y="Cell Count", title="Distribution of Somite Count Prediction Scores") +
  annotate("text", x=score_cutoff, y=Inf, label=paste0("20th percentile: ", round(score_cutoff, 3)), vjust=2, hjust=0, color="red")
dev.off()

# Day prediction score
score_cutoff <- quantile(df$predicted.day.score, 0.2, na.rm = TRUE)
pdf('images/scRNAseq/OMG/Distr_PredictionScore_day_20.pdf', width=8, height=5)
ggplot(df, aes(x=predicted.day.score)) +
  geom_histogram(bins=50, fill="#2ca02c", color="black") +
  geom_vline(xintercept = score_cutoff, color="red", linetype="dashed", size=1) +
  theme_classic() +
  labs(x="Predicted Day Score", y="Cell Count", title="Distribution of Day Prediction Scores") +
  annotate("text", x=score_cutoff, y=Inf, label=paste0("20th percentile: ", round(score_cutoff, 3)), vjust=2, hjust=0, color="red")
dev.off()

# apply prediction score filter 
percentile_celltype <- quantile(df$predicted.celltype_updated.score, 0.2, na.rm = TRUE)
percentile_somite <- quantile(df$predicted.somite_count.score, 0.2, na.rm = TRUE)
percentile_day <- quantile(df$predicted.day.score, 0.2, na.rm = TRUE)

Asmb <- subset(Asmb, subset = predicted.celltype_updated.score >= percentile_celltype & predicted.somite_count.score >= percentile_somite & predicted.day.score >= percentile_day)
TLS_projected <- subset(TLS_projected, subset = predicted.celltype_updated.score >= percentile_celltype & predicted.somite_count.score >= percentile_somite & predicted.day.score >= percentile_day)

combined <- merge(Asmb, TLS_projected)
df <- combined@meta.data
table(Asmb$condition, Asmb$predicted.celltype_updated)          

orig_df <- df 
ct10_by_cond <- orig_df %>%
  group_by(condition, predicted.celltype_updated) %>%
  tally(name = "n") %>%
  filter(n >= 10) %>%
  arrange(condition, desc(n))

print(ct10_by_cond, n=41)

allowed_pairs <- ct10_by_cond %>%
  mutate(key = paste(condition, predicted.celltype_updated, sep = "||")) %>%
  pull(key)

meta_asmb <- Asmb@meta.data
meta_asmb$key <- paste(meta_asmb$condition, meta_asmb$predicted.celltype_updated, sep = "||")
keep_asmb <- rownames(meta_asmb)[meta_asmb$key %in% allowed_pairs]
Asmb <- subset(Asmb, cells = keep_asmb)

meta_tls <- TLS_projected@meta.data
meta_tls$key <- paste(meta_tls$condition, meta_tls$predicted.celltype_updated, sep = "||")
keep_tls <- rownames(meta_tls)[meta_tls$key %in% allowed_pairs]
TLS_projected <- subset(TLS_OMG, cells = keep_tls)

table(Asmb$condition, Asmb$predicted.celltype_updated)
table(TLS_projected$condition, TLS_projected$predicted.celltype_updated)
     
# renormalize, rescale and rerun pca, umap 
Asmb <- NormalizeData(Asmb)
all.genes <- rownames(Asmb)
Asmb <- ScaleData(Asmb, features = all.genes)
Asmb <- RunPCA(Asmb, npcs = 50, features = VariableFeatures(Asmb))
Asmb <- RunUMAP(Asmb, dims = 1:30, return.model=TRUE) 

pdf("images/scRNAseq/Asmb_umap_OMG_celltype.pdf", width = 20, height = 6)
DimPlot(Asmb)+
FeaturePlot(Asmb, features = "percent.mt")+
DimPlot(Asmb, group.by = "predicted.celltype_updated",  cols = celltype_updated_colors)
dev.off()

pdf("images/scRNAseq/Asmb_umap_OMG_day.pdf", width = 14, height = 4)
DimPlot(Asmb)+
DimPlot(Asmb, group.by = "predicted.day", cols = day_colors)
dev.off()

# stacked barplot of cell type distribution across conditions
celltype_dist <- df %>%
  dplyr::group_by(condition, predicted.celltype_updated) %>%
  summarise(count = n()) %>%
  group_by(condition) %>%
  mutate(percent = count / sum(count) * 100)

pdf("images/scRNAseq/OMG/Asmb_celltype_distribution.pdf", width = 10, height = 6)
ggplot(celltype_dist, aes(x = condition, y = percent, fill = predicted.celltype_updated)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  # use color palette from celltype_updated_colors
  scale_fill_manual(values = celltype_updated_colors) +
  labs(x = "Condition", y = "Percentage of Cells", fill = "Predicted Cell Type", title = "Distribution of Predicted Cell Types Across Conditions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# plot marker expression 
# for each condition, calculate top markers per predicted.celltype_updated and plot heatmap doheatmap for each condition of top 5 markers per cell state
for(c in unique(Asmb$condition)){
  a <- subset(Asmb, condition == c)
  Idents(a) <- a$predicted.celltype_updated
  markers <- FindAllMarkers(a, only.pos=TRUE,   min.pct = 0.05, logfc.threshold = 0.25)
  markers <- subset(markers, p_val_adj < 0.05 & avg_log2FC > 1)
  top5 <- markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  pdf(paste0("images/scRNAseq/OMG/Asmb_marker_heatmap_", c, ".pdf"), width = 10, height = 6)
  print(DoHeatmap(a, features = top5$gene,  group.colors= celltype_updated_colors, size=3))
  dev.off()
}

known_markers <- c(
  "Foxg1", "Six3", "Lhx2", "Otx1", "Otx2", "Emx2", "Pax6", "Dmbx1", "Fgf8", "Wnt1",
  "En1", "En2", "Pax2", "Pax5", "Pax8", "Egr2", "Mafb", "Gbx2", "Hoxb1", "Pax3",
  "Rfx4", "Zic1", "Olig2", "Dbx1", "Hes5", "Dbx2", "Nkx6-1", "T", "Sox2", "Hes7"
)

for(c in unique(Asmb$condition)){
  a <- subset(Asmb, condition == c)
  Idents(a) <- a$predicted.celltype_updated
  pdf(paste0("images/scRNAseq/OMG/Asmb_known_marker_heatmap_", c, ".pdf"), width = 10, height = 6)
  print(DoHeatmap(a, features = known_markers,  group.colors= celltype_updated_colors, size=3))
  dev.off()
}

# FeaturePlot for Otx2
pdf("images/scRNAseq/OMG/Otx2_featureplot.pdf", width=14, height=6)
DimPlot(Asmb, group.by = "predicted.celltype_updated", cols = celltype_updated_colors)+
FeaturePlot(Asmb, features = "Otx2")
dev.off()

# save Asmb object 
saveRDS(Asmb, file = "data/scRNAseq/Asmb_OMG.rds")
saveRDS(TLS_projected, file = "data/scRNAseq/TLS_projected_OMG.rds")
Asmb <- readRDS("data/scRNAseq/Asmb_OMG.rds")

HAP <- subset(Asmb, condition == "Hypoxic")
HAP <- NormalizeData(HAP)
all.genes <- rownames(HAP)
HAP <- ScaleData(HAP, features=all.genes)
HAP <- RunPCA(HAP, npcs=50, features = VariableFeatures(HAP))
HAP <- RunUMAP(HAP, dims = 1:30, return.model=TRUE)

# save HAP object
saveRDS(HAP, file = "data/scRNAseq/HAP_OMG.rds")
HAP <- readRDS("data/scRNAseq/HAP_OMG.rds")

asmb_ct <- unique(Asmb$predicted.celltype_updated)
OMG_filtered <- subset(OMG, celltype_updated %in% asmb_ct)

# process OMG_filtered
OMG_filtered <- NormalizeData(OMG_filtered, verbose = FALSE)
OMG_filtered <- FindVariableFeatures(OMG_filtered, verbose = FALSE)
OMG_filtered <- ScaleData(OMG_filtered, verbose = FALSE)
OMG_filtered <- RunPCA(OMG_filtered, npcs = 50, verbose = FALSE)
OMG_filtered <- RunUMAP(OMG_filtered, dims = 1:30, return.model = T)

pdf('images/scRNAseq/OMG/UMAP_filtered_OMG_20f_celltype.pdf', width=12)
DimPlot(OMG_filtered, pt.size=1.5, reduction = "umap", group.by='celltype_updated', raster=F, cols=celltype_updated_colors)+coord_fixed(ratio=1)
dev.off()

pdf('images/scRNAseq/OMG/UMAP_filtered_OMG_day.pdf', width=8)
DimPlot(OMG_filtered, pt.size=1.5, reduction = "umap", group.by='day', raster=F, cols=day_colors, shuffle=TRUE)+coord_fixed(ratio=1)
dev.off()

## map Asmb onto OMG_filtered
Asmb_OMG <- Asmb
TLS_OMG <- TLS_projected

Asmb_OMG <- NormalizeData(Asmb_OMG, verbose = FALSE)
Asmb_OMG <- FindVariableFeatures(Asmb_OMG, verbose = FALSE)
features <- rownames(GetAssayData(Asmb_OMG, assay = "RNA", slot = "data"))
Asmb_OMG <- ScaleData(Asmb_OMG, features = features, verbose = FALSE)
Asmb_OMG <- RunPCA(Asmb_OMG, npcs = 50, verbose = FALSE)

TLS_OMG <- NormalizeData(TLS_OMG, verbose = FALSE)  
TLS_OMG <- FindVariableFeatures(TLS_OMG, verbose = FALSE)
all_features_tls <- rownames(GetAssayData(TLS_OMG, assay = "RNA", slot = "data"))
TLS_OMG <- ScaleData(TLS_OMG,  features = all_features_tls,verbose = FALSE)
TLS_OMG <- RunPCA(TLS_OMG, npcs = 50, verbose = FALSE)

anchors_Asmb <- FindTransferAnchors(reference = OMG_filtered, query = Asmb_OMG, dims = 1:50, reference.reduction = "pca")  
anchors_TLS_norm <- FindTransferAnchors(reference = OMG_filtered, query = TLS_OMG, dims = 1:50, reference.reduction = "pca")
Asmb_OMG <- MapQuery(anchorset = anchors_Asmb, reference = OMG_filtered, query = Asmb_OMG,  reference.reduction = "pca", reduction.model = "umap")
TLS_OMG <- MapQuery(anchorset = anchors_TLS_norm, reference = OMG_filtered, query = TLS_OMG, reference.reduction = "pca", reduction.model = "umap")

umap_data_new <- as.data.frame(Asmb_OMG[["ref.umap"]]@cell.embeddings)
head(umap_data_new)
colnames(umap_data_new) <- c("refumap_1", "refumap_2")
umap_data_WT <- as.data.frame(OMG_filtered[["umap"]]@cell.embeddings)
head(umap_data_WT) ##  umap_1 
umap_data_old <- as.data.frame(TLS_OMG[["ref.umap"]]@cell.embeddings)
head(umap_data_old)

umap_data_old$condition <- TLS_OMG@meta.data[row.names(umap_data_old),"condition"]
umap_data_new$condition <- Asmb_OMG@meta.data[row.names(umap_data_new),"condition"]
umap_data_old$cell_state_annot <- TLS_OMG@meta.data[row.names(umap_data_old),"predicted.celltype_updated"]
umap_data_new$cell_state_annot <- Asmb_OMG@meta.data[row.names(umap_data_new),"predicted.celltype_updated"]
umap_data_old$day <- TLS_OMG@meta.data[row.names(umap_data_old),"predicted.day"]
umap_data_new$day <- Asmb_OMG@meta.data[row.names(umap_data_new),"predicted.day"]

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
unique(umap_combined$cell_state_annot)

umap_combined$cell_state_annot <- factor(
  umap_combined$cell_state_annot,
  levels = names(cell_type_colored)
)
umap_combined$cell_type_num <- unname(
  celltype_order[ as.character(umap_combined$cell_state_annot) ]
)
umap_combined$cell_type_num <- factor(
  umap_combined$cell_type_num,
  levels = unname(celltype_order)   # "Neuromesodermal progenitors (1)", ...,"
)
head(umap_combined)
umap_combined$num <- str_extract(umap_combined$cell_type_num, "\\d+")


# Loop through each condition and save a separate plot
for (cond in conditions) {
  
  # Filter data for the current condition
  umap_subset <- subset(umap_combined, condition == cond)
  
  centroids <- umap_subset %>%
    group_by(cell_type_num, num) %>%
    summarise(refumap_1 = median(refumap_1), refumap_2 = median(refumap_2), .groups = "drop")

  # Generate plot
  single_plot <- ggplot() +
    rasterise(
      geom_point(data = umap_data_WT_rep, aes(x = umap_1, y = umap_2), 
                 size = 3, colour = "grey82", alpha = 0.5),
      dpi = 300
    ) +
    geom_point(data = umap_subset, aes(x = refumap_1, y = refumap_2, color = cell_type_num), 
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
    # add num (umap_combined$num) once to the plot for each cell type 
    geom_text(data = centroids, aes(x = refumap_1, y = refumap_2, label = num),
              size = 6, fontface = "bold", color = "black") +
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
  
  ggsave(paste0("images/scRNAseq/OMG/UMAP_", cond, "_projected_OMG.pdf"),
         single_plot, width = 10, height = 7)

  print(cond)
  print(dim(umap_subset))
}

# umap for wildtype reference 
umap_data_WT$cell_type <-  OMG_filtered@meta.data[row.names(umap_data_WT),"celltype_updated"]
umap_data_WT$cell_type_num <- unname(
  celltype_order[ as.character(umap_data_WT$cell_type) ]
)
umap_data_WT$cell_type_num <- factor(
  umap_data_WT$cell_type_num,
  levels = unname(celltype_order)   # "Neuromesodermal progenitors (1)", ..., "Unassigned 2 (20)"
)
umap_data_WT$num <- str_extract(umap_data_WT$cell_type_num, "\\d+")
# centroids for labels
centroids_WT <- umap_data_WT %>%
  group_by(cell_type_num, num) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p<- ggplot() +
    geom_point(data = umap_data_WT, aes(x = umap_1, y = umap_2, color = cell_type_num), 
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
    geom_text(data = centroids_WT, aes(x = umap_1, y = umap_2, label = num),
              size = 6, fontface = "bold", color = "black") +
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
  
ggsave(paste0("images/scRNAseq/OMG/UMAP_ref_OMG.pdf"), p, width = 10, height=7)

Asmb_OMG$condition <- ifelse(Asmb_OMG$condition == "Normoxic", "NAP", 
                             ifelse(Asmb_OMG$condition == "Hypo_XAV", "HAP+X", 
                                    ifelse(Asmb_OMG$condition == "Hypoxic", "HAP", "HAP Hif1a KO")))

calculate_cell_state_percentage <- function(seurat_object) {
  seurat_object@meta.data %>%
    group_by(condition, predicted.celltype_updated) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(condition) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    arrange(condition, desc(percentage))
}

# Compute cell state percentages for both objects
tls_plot_data <- calculate_cell_state_percentage(TLS_OMG)
asmb_plot_data <- calculate_cell_state_percentage(Asmb_OMG)
plot_data <- bind_rows(tls_plot_data, asmb_plot_data)

# Ensure ordered levels for plotting
all_levels     <- names(cell_type_colored)
present_levels <- intersect(all_levels, unique(plot_data$predicted.celltype_updated))
plot_data$predicted.celltype_updated <- factor(plot_data$predicted.celltype_updated,
                                              levels = present_levels)

# Plot
pdf("images/scRNAseq/OMG/stacked_OMG_predicted_states_conditions.pdf", width=14)
ggplot(plot_data, aes(x = factor(condition), y = percentage, fill = predicted.celltype_updated)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = cell_type_colored,
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

calculate_cell_state_percentage_binomCI <- function(seurat_object) {
  seurat_object@meta.data %>%
    group_by(condition, predicted.celltype_updated) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(condition) %>%
    mutate(n_cond = sum(count)) %>%
    rowwise() %>%
    mutate(percentage = count / n_cond * 100,
           ci_lo = binom.test(count, n_cond)$conf.int[1] * 100,
           ci_hi = binom.test(count, n_cond)$conf.int[2] * 100) %>%
    ungroup() %>%
    arrange(condition, desc(percentage))
}

tls_plot_data  <- calculate_cell_state_percentage_binomCI(TLS_OMG)
asmb_plot_data <- calculate_cell_state_percentage_binomCI(Asmb_OMG)
plot_data <- bind_rows(tls_plot_data, asmb_plot_data)

all_levels     <- names(cell_type_colored)
present_levels <- intersect(all_levels, unique(plot_data$predicted.celltype_updated))
plot_data$predicted.celltype_updated <- factor(plot_data$predicted.celltype_updated,
                                                levels = present_levels)

pdf("images/scRNAseq/OMG/OMG_predicted_states_conditions_CI.pdf",
    width = 16, height = 10)
ggplot(plot_data,
       aes(x = factor(condition), y = percentage,
           fill = predicted.celltype_updated)) +
  geom_col() +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.3) +
  facet_wrap(~ predicted.celltype_updated, scales = "free_y") +
  scale_fill_manual(values = cell_type_colored, breaks = present_levels) +
  labs(x = "", y = "Percentage of cells (95% binomial CI)", fill = "Cell states") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(size = 11, angle = 45, hjust = 1, color = "black"),
    axis.text.y  = element_text(size = 10, color = "black"),
    strip.background = element_blank(),
    strip.text   = element_text(size = 9),
    legend.position = "none"   # facet labels name the types; legend redundant
  )
dev.off()

## Staging, dotplot for HAPs
combined <- merge(Asmb_OMG, TLS_OMG)
dim(combined)
df <- combined@meta.data
haps <- subset(df, condition == "HAP")

# DotPlot to show cell fraction per predicted stage 
dotplot_data <- haps %>%
  dplyr::group_by(condition, predicted.day) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%  # normalization within each condition
  dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) 
print(dotplot_data)

all_stages <- c("E8.0-E8.5","E8.75", "E9.0","E9.25","E9.5","E9.75")
# ensure order of stages in plot 
dotplot_data$predicted.day <- factor(dotplot_data$predicted.day, levels = all_stages, ordered = TRUE)

cond <- "HAP"
# dot plots for each condition pair
for (pair in cond) {
  
  dotplot_subset <- haps %>%
    filter(cond %in% pair) %>%
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
      legend.title = element_text(size = 14), aspect.ratio = 0.5
    )
  filename <- paste0("images/scRNAseq/OMG/OMG_dotplot_stage_", pair[1], ".pdf")
  ggsave(filename, width = 7, height = 2.8, plot)
  print(pair)
}

# turn somite count into numeric number 
haps$predicted.somite_num <- as.numeric(gsub(" somites", "", haps$predicted.somite_count))
all_somites <- min(haps$predicted.somite_num, na.rm=TRUE) : max(haps$predicted.somite_num, na.rm=TRUE)

# dot plots for each condition pair
for (pair in cond) {
  
  dotplot_subset <- haps %>%
    filter(cond %in% pair) %>%
    dplyr::group_by(condition, predicted.somite_num) %>%
    dplyr::summarise(cell_count = n(), .groups = "drop") %>%
    dplyr::group_by(condition) %>%  
    dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) %>%
    ungroup()  
  dotplot_subset <- dotplot_subset %>%
    tidyr::complete(condition = pair, predicted.somite_num = all_somites, fill = list(cell_count = NA, percentage = NA))
  dotplot_subset$predicted.somite_num <- factor(dotplot_subset$predicted.somite_num, levels = all_somites, ordered = TRUE)
  
  plot <- ggplot(dotplot_subset, aes(x = predicted.somite_num, y = condition, size = percentage)) +
    geom_point(alpha = 0.8, na.rm = TRUE) + 
    scale_size_continuous(range = c(2, 10)) +  # No dots for 0% values
    scale_x_discrete(expand = c(0.02, 0.02)) +
     coord_cartesian(clip = "off") +
    labs(x = "Predicted somite count", y = "", size = "Cell fraction (%)", title = "") +
    theme_minimal() +
    theme_classic() +
    theme(
      # rotate x-axis labels 
      axis.text.x = element_text(size = 14, color = "black", angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),  
      axis.title = element_text(size = 16),  
      legend.text = element_text(size = 12),  
      legend.title = element_text(size = 14), aspect.ratio = 0.1
    )
  filename <- paste0("images/scRNAseq/OMG/OMG_dotplot_somites_", pair[1], ".pdf")
  ggsave(filename, width = 11, height = 2, plot)
  print(pair)
}

somite_scores <- GetAssayData(Asmb_OMG, assay = "prediction.score.somite_count", slot = "data")
hapa_cells <- rownames(haps)
somite_scores <- somite_scores[, hapa_cells]

# get first and seconds top predictions for heatmap
top1 <- apply(somite_scores, 2, function(x) {
  names(which.max(x))
})
top2 <- apply(somite_scores, 2, function(x) {
  names(sort(x, decreasing = TRUE)[2])
})

# plot heatmap of top1 vs top2 predictions 
df <- data.frame(
  top1 = as.character(top1),
  top2 = as.character(top2)
)
df$top1_num <- as.numeric(gsub(" somites", "", df$top1))
df$top2_num <- as.numeric(gsub(" somites", "", df$top2))
head(df)

lvl <- 0:31
tbl <- with(df, table(factor(top1_num, levels = lvl),
                      factor(top2_num, levels = lvl)))
mat <- as.matrix(tbl)
# Normalize to 0–1 by max 
den <- max(mat)
mat_norm <- if (den > 0) mat/den else mat
lvl <- as.character(0:31)
mat_plot <- mat_norm[ lvl, lvl, drop = FALSE ]          
    
col_fun <- colorRamp2(c(0,  1), c("oldlace", "darkred"))

ht <- Heatmap(
  mat_plot,
  name = "Cell rate",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = lvl,                
  column_order = lvl,              
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_title = "Most likely predicted somite count",
  column_title = "Second most likely predicted somite count",
  heatmap_legend_param = list(title = "Cell fraction"),
  row_names_side = "left",
  column_title_side = "bottom",
  column_names_rot = 90
)
pdf("images/scRNAseq/OMG/OMG_heatmap_somite_top1_top2_HAPs.pdf",
    width = 7, height = 6)
draw(ht)
dev.off()

# Weighted sum (expected somite count) per cell
w <- as.numeric(gsub(" somites", "", rownames(somite_scores)))
num   <- colSums(sweep(somite_scores, 1, w, `*`))      # sum (score * somite)
colSums(somite_scores)                       
expected_somite <- num 

df_w <- data.frame(
  cell = colnames(somite_scores),
  expected_somite = expected_somite,
  row.names = colnames(somite_scores)
)

# Single boxplot (all cells)
pdf("images/scRNAseq/OMG/OMG_boxplot_somite_weighted_sum_HAPs.pdf",
       width = 2, height = 5)
ggplot(df_w, aes(y = expected_somite)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(y = "Predicted somite count (HAPs)") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16)
  )
dev.off()

df_sum <- data.frame(
  group = "HAP",
  mean  = mean(df_w$expected_somite, na.rm = TRUE),
  sd    = sd(df_w$expected_somite,   na.rm = TRUE)
)

# bar plot 
pdf("images/scRNAseq/OMG/OMG_barplot_somite_weighted_sum_HAPs.pdf",
    width = 2.5, height = 5)
  ggplot() +
  # bar at mean
  geom_col(data = df_sum, aes(x = group, y = mean),
           fill = "#963752", width = 0.6) +
  # SD error bar
  geom_errorbar(data = df_sum,
                aes(x = group, ymin = mean - sd, ymax = mean + sd),
                width = 0.2) +
  geom_jitter(data = df_w, aes(x = "HAP", y = expected_somite),
              width = 0.15, size = 0.5, alpha = 0.3, color = "black") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(y = "Predicted somite count (HAPs)") +
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16)
  )
dev.off()


# Staging per cell type 
day_scores <- GetAssayData(Asmb_OMG, assay = "prediction.score.day", slot = "data")
hapa_cells <- rownames(haps)
day_scores <- day_scores[, hapa_cells]

w <- sapply(rownames(day_scores), function(x) {
  nums <- as.numeric(regmatches(x, gregexpr("\\d+\\.?\\d*", x))[[1]])
  if (length(nums) == 2) mean(nums) else nums
})

# Weighted sum (expected somite count) per cell
num   <- colSums(sweep(day_scores, 1, w, `*`))      # sum (score * somite)
denom <- colSums(day_scores)                      
expected_day <- num / ifelse(denom == 0, 1, denom)

df_w <- data.frame(
  cell = colnames(day_scores),
  expected_day = expected_day,
  row.names = colnames(day_scores)
)

# add cell type
df_w$cell_type <- haps[colnames(day_scores), "predicted.celltype_updated"]
df_w$cell_type <- factor(df_w$cell_type,
                         levels = names(cell_type_colored))

yticks <- sort(unique(c(as.numeric(w), 8.5)))  

pdf("images/scRNAseq/OMG/OMG_boxplot_day_weighted_sum_HAPs_per_celltype.pdf",
    width = 12, height = 6)
ggplot(df_w, aes(x = cell_type, y = expected_day,  fill = cell_type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(y = "Predicted stage", x = "") +
   scale_fill_manual(values = cell_type_colored) +
   scale_y_continuous(
    breaks = yticks
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()

neural_ct <- c(
  "Anterior floor plate",
  "Floorplate and p3 domain",
  "Anterior roof plate",
  "Posterior roof plate",
  "Spinal cord ventral progenitors",
  "Spinal cord motor neurons",
  "Spinal cord dorsal progenitors",
  "Spinal cord/r7/r8",
  "NMPs and spinal cord progenitors",
  "Dorsal telencephalon",
  "Telencephalon",
  "Posterior Forebrain / Diencephalon",
  "Midbrain",
  "Hindbrain",
  "Midbrain-hindbrain boundary",
  "Hypothalamus",
  "Hypothalamus (Sim1+)",
  "Eye field",
  "Cranial motor neurons",
  "Neural crest (PNS neurons)",
  "Neural crest (PNS glia)",
  "Glutamatergic neurons",
  "Neural progenitor cells (Neurod1+)"
)

Asmb_neural <- subset(Asmb_OMG, subset = predicted.celltype_updated %in% neural_ct)
neural_asmb <- unique(Asmb_neural$predicted.celltype_updated)
OMG_neural <- subset(OMG_filtered, subset = celltype_updated %in% neural_asmb)
unique(OMG_neural$celltype_updated)
# [1] "Floorplate and p3 domain"           "Midbrain-hindbrain boundary"       
# [3] "NMPs and spinal cord progenitors"   "Spinal cord/r7/r8"                 
# [5] "Hindbrain"                          "Posterior Forebrain / Diencephalon"
# [7] "Telencephalon"                      "Midbrain"  


umap_data_new <- as.data.frame(Asmb_OMG[["ref.umap"]]@cell.embeddings)
colnames(umap_data_new) <- c("refumap_1", "refumap_2")
head(umap_data_new)
umap_data_WT <- as.data.frame(OMG_filtered[["umap"]]@cell.embeddings)
head(umap_data_WT) ##  umap_1 

umap_data_new$condition <- Asmb_OMG@meta.data[row.names(umap_data_new),"condition"]
umap_data_new$cell_state_annot <- Asmb_OMG@meta.data[row.names(umap_data_new),"predicted.celltype_updated"]
umap_data_new$day <- Asmb_OMG@meta.data[row.names(umap_data_new),"predicted.day"]
umap_combined <- bind_rows(
  umap_data_new %>% mutate(source = "Asmb")
)
umap_combined$condition <- factor(umap_combined$condition)

umap_data_WT_rep <- umap_data_WT[rep(1:nrow(umap_data_WT), times = length(unique(umap_combined$condition))), ]
umap_data_WT_rep$condition <- rep(unique(umap_combined$condition), each = nrow(umap_data_WT))

label_positions <- umap_combined %>%
  group_by(cell_state_annot) %>%
  summarise(refumap_1 = mean(refumap_1), refumap_2 = mean(refumap_2)) %>%
  mutate(cell_state_id = str_extract(cell_state_annot, "\\d+"))  # Extract only numbers

conditions <- unique(umap_combined$condition)

umap_combined$cell_state_annot <- factor(
  umap_combined$cell_state_annot,
  levels = names(cell_type_colored)
)

umap_data_WT$day <- OMG_filtered@meta.data[row.names(umap_data_WT),"day"]
all_days <- c("E8.0-E8.5","E8.75", "E9.0","E9.25","E9.5","E9.75")
umap_data_WT$day <- factor(umap_data_WT$day, levels = all_days, ordered = TRUE)
umap_data_WT <- umap_data_WT[sample(nrow(umap_data_WT)), ]

p <- ggplot() +
    geom_point(data = umap_data_WT, aes(x = umap_1, y = umap_2, color = day), 
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
    guides(color = guide_legend(title = "", override.aes = list(size = 5))) +
    coord_fixed(ratio = 1) 
ggsave(paste0("images/scRNAseq/OMG/UMAP_ref_day_OMG.pdf"), p, width = 10, height=7)


table(Asmb_neural$condition)

# Compute cell state percentages for both objects
asmb_plot_data <- calculate_cell_state_percentage(Asmb_neural)
plot_data <- subset(asmb_plot_data, condition %in% c("HAP", "HAP+X"))
all_levels     <- names(cell_type_colored)
present_levels <- intersect(all_levels, unique(plot_data$predicted.celltype_updated))
plot_data$predicted.celltype_updated <- factor(plot_data$predicted.celltype_updated,
                                              levels = present_levels)

pdf("images/scRNAseq/OMG/stacked_OMG_predicted_states_neural_85.pdf", width=9)
ggplot(plot_data, aes(x = factor(condition), y = percentage, fill = predicted.celltype_updated)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colored) + 
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
    legend.title = element_text(size = 16, color = "black"))
dev.off()


neural_main <- c("Posterior Forebrain / Diencephalon", "Midbrain", "Midbrain-hindbrain boundary", "Hindbrain", "Spinal cord/r7/r8", "NMPs and spinal cord progenitors")

# Define condition pairs for comparisons
condition_pairs <- list(
  c("HAP", "HAP+X")
)

# marker groups
marker_groups <- list(
  "Posterior Forebrain / Diencephalon" = c("Foxg1", "Six3", "Lhx2", "Otx1", "Emx2", "Pax6", "Fezf1", "Rax"),
  "Midbrain" = c("Dmbx1", "Otx2"),
  "Midbrain-hindbrain boundary" = c("Fgf8", "Wnt1", "En1", "En2", "Pax2", "Pax5", "Pax8"),
  "Hindbrain" = c("Egr2", "Mafb", "Hoxb1", "Hoxa2", "Gbx2"),
  "Spinal cord/r7/r8" = c("Zic1", "Olig2", "Dbx1", "Hes5", "Dbx2", "Rfx4", "Nkx6-1", "Nkx6-2", "Pax3"),
  "NMPs and spinal cord progenitors" = c("T", "Sox2", "Hes")
)

markers <- unlist(marker_groups, use.names = FALSE)

Asmb_neural <- subset(Asmb_neural, subset = predicted.celltype_updated %in% neural_main)
Asmb_neural <- ScaleData(Asmb_neural, features= markers, verbose = FALSE)

create_dotplot <- function(seurat_obj, condition_pair, markers) {
  obj <- subset(seurat_obj, condition %in% condition_pair)
  obj <- subset(obj, predicted.celltype_updated %in% neural_main)

  # Run DotPlot per condition 
  df <- do.call(rbind, lapply(condition_pair, function(cond) {
    sub_obj <- subset(obj, condition == cond)
    Idents(sub_obj) <- "predicted.celltype_updated"
    d <- DotPlot(sub_obj,
                 features  = markers,
                 group.by  = "predicted.celltype_updated",
                 dot.scale = 8)$data
    d$condition <- cond
    d$celltype  <- as.character(d$id)
    d
  }))

  df$condition <- factor(df$condition, levels = condition_pair)

  ggplot(df, aes(x = features.plot, y = celltype)) +
    geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
    scale_size(range = c(0, 8), name = "% expressed") +
    scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red",
                           midpoint = 0,
                           name = "Avg. expression\n(scaled)") +
    facet_wrap(~ condition, ncol = 1) +
    theme_classic()+
     theme(axis.text.x      = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL)
}

for (pair in condition_pairs) {
  filename <- paste0("images/scRNAseq/OMG/marker_expression_",pair[1], "_vs_", pair[2], ".pdf")
  pdf(file = filename, width = 14, height = 8)
  p <- create_dotplot(Asmb_neural, pair, markers)
  print(p)
  dev.off()
  message("Saved plot for ", pair[1], " vs ", pair[2], " as ", filename)
}


test_pair <- c("HAP+X", "HAP")

subset <- subset(Asmb_neural, condition == test_pair[1] | condition == test_pair[2])
subset <- subset(subset, predicted.celltype_updated %in% neural_main)
subset <- ScaleData(subset, features= markers, verbose = FALSE)
subset$predicted.celltype_updated <- factor(
  subset$predicted.celltype_updated,
  levels = c("Posterior Forebrain / Diencephalon", "Midbrain", "Midbrain-hindbrain boundary", "Hindbrain",  "Spinal cord/r7/r8", "NMPs and spinal cord progenitors")
)
Idents(subset) <- "predicted.celltype_updated"

all_markers <- unique(unlist(marker_groups, use.names = FALSE))  # or just `markers`

scale_mat <- subset@assays$RNA@scale.data
scale_mat <- GetAssayData(subset, assay = "RNA", layer = "scale.data")
present <- intersect(all_markers, rownames(scale_mat))
missing <- setdiff(all_markers, rownames(scale_mat))

# zero matrix for missing genes
if (length(missing) > 0) {
  zero_mat <- matrix(
    0,
    nrow = length(missing),
    ncol = ncol(scale_mat),
    dimnames = list(missing, colnames(scale_mat))
  )
  subset_data <- rbind(scale_mat[present, , drop = FALSE], zero_mat)
} else {
  subset_data <- scale_mat[present, , drop = FALSE]
}

gene_order <- unique(unlist(marker_groups, use.names = FALSE))
subset_data <- subset_data[gene_order[gene_order %in% rownames(subset_data)], , drop = FALSE]

subset_meta <- subset@meta.data
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell") 

# Convert to long format
df <- as.data.frame(t(as.matrix(subset_data)))
df$cell <- rownames(df)
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell")  
df <- df %>%
  pivot_longer(cols = -cell, names_to = "gene", values_to = "expression") %>%
  left_join(subset_meta, by = "cell") 

# mean expression per group
df_summary <- df %>%
  group_by(gene,  condition, predicted.celltype_updated) %>%
  summarize(
    avg_expression = mean(expression),  
    num_expressing = sum(expression > 0),  # Count of cells expressing the gene
    total_cells = n(),  # Total cells in this group
    pct_expressed = (num_expressing / total_cells) * 100,  
    .groups = "drop"
  )

subset <- df_summary
gene_order <- unlist(marker_groups, use.names = FALSE)  
subset$gene <- factor(subset$gene, levels = gene_order)

filename <- paste0("images/scRNAseq/OMG/marker_expression_",test_pair[1], "_vs_", test_pair[2], ".pdf")
  pdf(file = filename,height=8, width=11)
ggplot(subset, aes(x = gene, y = condition, size = pct_expressed)) +
  geom_point(aes(fill = avg_expression), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "gray93", high = "#b2182b", midpoint = 0) +
  scale_size(range = c(1, 8)) +
  facet_wrap(~predicted.celltype_updated, ncol = 1) +
  theme_minimal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    text = element_text(color = "black"),             
    axis.text = element_text(color = "black"),        
    axis.title = element_text(color = "black"),       
    strip.text = element_text(color = "black"),       
    legend.text = element_text(color = "black"),      
    legend.title = element_text(color = "black")      
  ) +
  labs(x = "", y = "", fill = "Avg. Expression", size = "Percent Expressed") 
dev.off()


# Two-sample proportion test: neural cell-state composition, HAP+X vs HAP 
# Counts per state per condition (within the neural subset already built above)
state_counts <- subset@meta.data %>%   # `subset` = the two-condition neural object
  dplyr::count(condition, predicted.celltype_updated, name = "count") %>%
  tidyr::complete(condition, predicted.celltype_updated, fill = list(count = 0))

# Per-condition totals (denominator = all neural cells in that condition)
cond_totals <- state_counts %>%
  group_by(condition) %>% summarise(total = sum(count), .groups = "drop")

# For each state: 2x2 of (in-state vs not) x (HAP+X vs HAP) -> prop.test
prop_results <- state_counts %>%
  left_join(cond_totals, by = "condition") %>%
  group_by(predicted.celltype_updated) %>%
  summarise(
    k_A   = count[condition == test_pair[1]],          # HAP+X in-state
    n_A   = total[condition == test_pair[1]],          # HAP+X total neural
    k_B   = count[condition == test_pair[2]],          # HAP in-state
    n_B   = total[condition == test_pair[2]],          # HAP total neural
    pct_A = k_A / n_A * 100,
    pct_B = k_B / n_B * 100,
    p_value = prop.test(c(k_A, k_B), c(n_A, n_B))$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

print(prop_results)


# Proportion bar plot with significance
plot_props <- prop_results %>%
  select(predicted.celltype_updated, pct_A, pct_B) %>%
  tidyr::pivot_longer(c(pct_A, pct_B),
                      names_to = "which", values_to = "percentage") %>%
  mutate(condition = ifelse(which == "pct_A", test_pair[1], test_pair[2]),
         condition = factor(condition, levels = test_pair))

state_levels <- c("Posterior Forebrain / Diencephalon", "Midbrain",
                  "Midbrain-hindbrain boundary", "Hindbrain",
                  "Spinal cord/r7/r8", "NMPs and spinal cord progenitors")
plot_props$predicted.celltype_updated <-
  factor(plot_props$predicted.celltype_updated, levels = state_levels)

sig <- prop_results %>%
  mutate(predicted.celltype_updated =
           factor(predicted.celltype_updated, levels = state_levels),
         label = ifelse(p_adj < 0.001, "***",
                 ifelse(p_adj < 0.01,  "**",
                 ifelse(p_adj < 0.05,  "*", ""))),
         y = pmax(pct_A, pct_B) + 2)

filename <- paste0("images/scRNAseq/OMG/neural_proportions_",
                   test_pair[1], "_vs_", test_pair[2], ".pdf")
pdf(filename, height = 6, width = 13)
ggplot(plot_props, aes(x = condition, y = percentage,
                       fill = predicted.celltype_updated)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(data = sig,
            aes(x = 1.5, y = y, label = label),      
            inherit.aes = FALSE, size = 5) +
  facet_wrap(~ predicted.celltype_updated, nrow = 1) +
  scale_fill_manual(values = celltype_updated_colors, guide = "none") +
  labs(x = "", y = "Percentage of neural cells") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9))
dev.off()

# only for condition == HAPA
test_pair <- c("HAP")
subset <- subset(Asmb_neural, condition == test_pair[1])
subset <- subset(subset, predicted.celltype_updated %in% neural_main)
subset <- ScaleData(subset, features= markers, verbose = FALSE)
subset$predicted.celltype_updated <- factor(
  subset$predicted.celltype_updated,
  levels = c("Posterior Forebrain / Diencephalon", "Midbrain", "Midbrain-hindbrain boundary", "Hindbrain",  "Spinal cord/r7/r8", "NMPs and spinal cord progenitors")
)
Idents(subset) <- "predicted.celltype_updated"

scale_mat <- GetAssayData(subset, assay = "RNA", layer = "scale.data")
present <- intersect(all_markers, rownames(scale_mat))
missing <- setdiff(all_markers, rownames(scale_mat))

# zero matrix for missing genes
if (length(missing) > 0) {
  zero_mat <- matrix(
    0,
    nrow = length(missing),
    ncol = ncol(scale_mat),
    dimnames = list(missing, colnames(scale_mat))
  )
  subset_data <- rbind(scale_mat[present, , drop = FALSE], zero_mat)
} else {
  subset_data <- scale_mat[present, , drop = FALSE]
}

gene_order <- unique(unlist(marker_groups, use.names = FALSE))
subset_data <- subset_data[gene_order[gene_order %in% rownames(subset_data)], , drop = FALSE]
subset_meta <- subset@meta.data
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell") 

df <- as.data.frame(t(as.matrix(subset_data)))
df$cell <- rownames(df)
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell")  
df <- df %>%
  pivot_longer(cols = -cell, names_to = "gene", values_to = "expression") %>%
  left_join(subset_meta, by = "cell") 

# mean expression per group
df_summary <- df %>%
  group_by(gene,  condition, predicted.celltype_updated) %>%
  summarize(
    avg_expression = mean(expression),  
    num_expressing = sum(expression > 0),  
    total_cells = n(),  
    pct_expressed = (num_expressing / total_cells) * 100,  
    .groups = "drop"
  )

subset <- df_summary
gene_order <- unlist(marker_groups, use.names = FALSE)  # Get all genes in order
subset$gene <- factor(subset$gene, levels = gene_order)

filename <- paste0("images/scRNAseq/OMG/marker_expression_",test_pair[1], ".pdf")
pdf(file = filename,height=4, width=14)
ggplot(subset, aes(x = gene, y = predicted.celltype_updated, size = pct_expressed)) +
  geom_point(aes(fill = avg_expression), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "gray93", high = "#b2182b", midpoint = 0) +
  scale_size(range = c(1, 8)) +
  scale_y_discrete(limits = rev) + 
  theme_minimal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    text = element_text(color = "black"),              # Sets all text to black
    axis.text = element_text(color = "black"),         # Black axis text
    axis.title = element_text(color = "black"),        # Black axis labels
    strip.text = element_text(color = "black"),        # Black facet labels
    legend.text = element_text(color = "black"),       # Black legend text
    legend.title = element_text(color = "black")       # Black legend title
  ) +
  labs(x = "", y = "", fill = "Avg. Expression", size = "Percent Expressed") 
dev.off()

all.genes <- rownames(Asmb_neural)
Asmb_neural <- ScaleData(Asmb_neural, features = all.genes)
mb_markers_omg <- FindMarkers(OMG_neural, ident.1 = "Midbrain", only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
mb_markers_omg <- subset(mb_markers_omg, p_val_adj < 0.05 & avg_log2FC > 1)
mb_markers_omg[rownames(mb_markers_omg) %in% c("Dmbx1", "Otx2"), ]
top20 <- rownames(head(mb_markers_omg[order(mb_markers_omg$avg_log2FC, decreasing = TRUE), ], 20))

test_pair <- c("HAP")
subset <- subset(Asmb_neural, condition == test_pair[1])
subset <- subset(subset, predicted.celltype_updated %in% neural_main)
subset <- ScaleData(subset, features= top20, verbose = FALSE)
subset$predicted.celltype_updated <- factor(
  subset$predicted.celltype_updated,
  levels = c("Posterior Forebrain / Diencephalon", "Midbrain", "Midbrain-hindbrain boundary", "Hindbrain",  "Spinal cord/r7/r8", "NMPs and spinal cord progenitors")
)
Idents(subset) <- "predicted.celltype_updated"

scale_mat <- GetAssayData(subset, assay = "RNA", layer = "scale.data")
present <- intersect(top20, rownames(scale_mat))
missing <- setdiff(top20, rownames(scale_mat))

# zero matrix for missing genes
if (length(missing) > 0) {
  zero_mat <- matrix(
    0,
    nrow = length(missing),
    ncol = ncol(scale_mat),
    dimnames = list(missing, colnames(scale_mat))
  )
  subset_data <- rbind(scale_mat[present, , drop = FALSE], zero_mat)
} else {
  subset_data <- scale_mat[present, , drop = FALSE]
}

subset_meta <- subset@meta.data
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell") 

df <- as.data.frame(t(as.matrix(subset_data)))
df$cell <- rownames(df)
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell")  
df <- df %>%
  pivot_longer(cols = -cell, names_to = "gene", values_to = "expression") %>%
  left_join(subset_meta, by = "cell") 

# mean expression per group
df_summary <- df %>%
  group_by(gene,  condition, predicted.celltype_updated) %>%
  summarize(
    avg_expression = mean(expression),  
    num_expressing = sum(expression > 0), 
    total_cells = n(),  
    pct_expressed = (num_expressing / total_cells) * 100,  
  )

subset <- df_summary
filename <- paste0("images/scRNAseq/OMG/top20_Midbrain_OMG_marker_expression",test_pair[1], ".pdf")
pdf(file = filename,height=5, width=14)
ggplot(subset, aes(x = gene, y = predicted.celltype_updated, size = pct_expressed)) +
  geom_point(aes(fill = avg_expression), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "gray93", high = "#b2182b", midpoint = 0) +
  scale_size(range = c(1, 8)) +
  scale_y_discrete(limits = rev) + 
  theme_minimal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    text = element_text(color = "black"),              
    axis.text = element_text(color = "black"),         
    axis.title = element_text(color = "black"),        
    strip.text = element_text(color = "black"),        
    legend.text = element_text(color = "black"),       
    legend.title = element_text(color = "black")       
  ) +
  labs(x = "", y = "", fill = "Avg. Expression", size = "Percent Expressed") 
dev.off()


Asmb_neural <- ScaleData(Asmb_neural, features= markers, verbose = FALSE)
# plot markers for OMG_filtered 
Idents(OMG_filtered) <- "celltype_updated"
markers <- FindAllMarkers(OMG_filtered, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
# filter markers
markers <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 1, ]
# write markers to csv & excel
write.csv(markers, "tables/scRNAseq/OMG/OMG_filtered_markers.csv")
write.xlsx(markers, "tables/scRNAseq/OMG/OMG_filtered_markers.xlsx")
markers_tome <- read.csv("tables/scRNAseq/TOME/TOME_filtered_markers.csv")
markers$reference <- "OMG"
markers_tome$reference <- "TOME" 
markers_combined <- rbind(markers, markers_tome)
write.xlsx(markers_combined, "tables/scRNAseq/OMG/OMG_TOME_filtered_markers_combined.xlsx")


top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DefaultAssay(OMG_filtered) <- "RNA"
# z-scores 
OMG_filtered <- ScaleData(
  OMG_filtered,
  assay   = "RNA",
  features = rownames(OMG_filtered[["RNA"]]),
  verbose = FALSE
)

top10 <- top10 %>%
  mutate(cluster = factor(cluster, levels = names(celltype_order))) %>%
  arrange(cluster, desc(avg_log2FC))
feature_list <- split(top10$gene, top10$cluster)
Idents(OMG_filtered) <- factor(
  Idents(OMG_filtered),
  levels = names(celltype_order) 
)

#heatmap
pdf("images/scRNAseq/OMG/Heatmap_marker_expression_OMG_filtered_top10.pdf", width=26, height=18)  
DoHeatmap(OMG_filtered, features = top10$gene, size=4, group.colors=cell_type_colored, raster=TRUE) + 
  theme(axis.text.y = element_text(size = 6))
dev.off()

# save all seurat objects
saveRDS(OMG_filtered, "data/scRNAseq/OMG_filtered.rds")
saveRDS(OMG_neural,"data/scRNAseq/OMG_neural.rds" )
saveRDS(Asmb_OMG,"data/scRNAseq/Asmb_OMG_filtered.rds"  )
saveRDS(TLS_OMG, "data/scRNAseq/TLS_OMG.rds")
saveRDS(Asmb_neural, "data/scRNAseq/Asmb_neural.rds" )
#########################################################################################################
# Gut markers 
# marker groups
marker_groups <- list(
  "Gut" = c("Sox17", "Foxa2", "Hhex", "Pdx1", "Gata4", "Gata6", "Cdx1", "Cdx2", "Cdx4", "Evx1")
)

markers <- unlist(marker_groups, use.names = FALSE)
Asmb_gut <- subset(Asmb_OMG, subset = predicted.celltype_updated %in% c("Gut"))

create_dotplot <- function(seurat_obj, condition_pair, markers) {
  # Subset conditions
  subset <- subset(seurat_obj, condition == condition_pair[1] | condition == condition_pair[2] )
  subset <- subset(subset, predicted.celltype_updated %in% neural_main)
  # dot plot
  Idents(subset) <- "predicted.celltype_updated"
  DotPlot(subset, features = markers, group.by = "predicted.celltype_updated", split.by ="condition", cols = c("blue", "red"), dot.scale = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
}

test_pair <- c("HAP")
subset <- subset(Asmb_gut, condition == test_pair[1])
subset <- subset(subset, predicted.celltype_updated %in% c("Gut"))
subset <- ScaleData(subset, features= markers, verbose = FALSE)
Idents(subset) <- "predicted.celltype_updated"

scale_mat <- GetAssayData(subset, assay = "RNA", layer = "scale.data")
present <- intersect(markers, rownames(scale_mat))
missing <- setdiff(markers, rownames(scale_mat))

# zero matrix for missing genes
if (length(missing) > 0) {
  zero_mat <- matrix(
    0,
    nrow = length(missing),
    ncol = ncol(scale_mat),
    dimnames = list(missing, colnames(scale_mat))
  )
  subset_data <- rbind(scale_mat[present, , drop = FALSE], zero_mat)
} else {
  subset_data <- scale_mat[present, , drop = FALSE]
}

gene_order <- unique(unlist(marker_groups, use.names = FALSE))
subset_data <- subset_data[gene_order[gene_order %in% rownames(subset_data)], , drop = FALSE]
subset_meta <- subset@meta.data
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell") 
df <- as.data.frame(t(as.matrix(subset_data)))
df$cell <- rownames(df)
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell")  
df <- df %>%
  pivot_longer(cols = -cell, names_to = "gene", values_to = "expression") %>%
  left_join(subset_meta, by = "cell") 

# mean expression per group
df_summary <- df %>%
  group_by(gene,  condition, predicted.celltype_updated) %>%
  summarize(
    avg_expression = mean(expression),  
    num_expressing = sum(expression > 0),  
    total_cells = n(),  
    pct_expressed = (num_expressing / total_cells) * 100, 
    .groups = "drop"
  )

subset <- df_summary
gene_order <- unlist(marker_groups, use.names = FALSE)  # Get all genes in order
subset$gene <- factor(subset$gene, levels = gene_order)

# DotPlot
filename <- paste0("images/scRNAseq/OMG/marker_expression_Gut_",test_pair[1], ".pdf")
pdf(file = filename,height=5, width=6)
ggplot(subset, aes(x = gene, y = predicted.celltype_updated, size = pct_expressed)) +
  geom_point(shape = 21,  color = "black", stroke = 0.5, aes(fill = avg_expression)) +
  scale_fill_gradient2(low = "#2166ac", mid = "gray93", high = "#b2182b", midpoint = 0,  labels = function(x) {
    txt <- formatC(x, format = "e", digits = 1)
    parse(text = gsub("e\\+?", " %*% 10^", txt))
  }) +
  scale_size(range = c(1, 8)) +
  theme_minimal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    text = element_text(color = "black"),              # Sets all text to black
    axis.text = element_text(color = "black"),         # Black axis text
    axis.title = element_text(color = "black"),        # Black axis labels
    strip.text = element_text(color = "black"),        # Black facet labels
    legend.text = element_text(color = "black"),       # Black legend text
    legend.title = element_text(color = "black")       # Black legend title
  ) +
  labs(x = "", y = "", fill = "Avg. Expression", size = "Percent Expressed") + coord_fixed(ratio = 9/7) 
dev.off()

###############################################################
# more marker genes for somites 
# marker groups
marker_groups <- list(
  "Somites" = c("Tbx6", "Hes7", "T", "Foxc1", "Foxc2", "Meox1", "Mesp2", "Uncx", "EphA4", "Tbx18", "EphrinB2", "Pax3", "Pax7", "Sim1", "Pax1", "Pax9", "Sox9")
)

markers <- unlist(marker_groups, use.names = FALSE)
Asmb_somites <- subset(Asmb_OMG, subset = predicted.celltype_updated %in% c("Mesodermal progenitors (Tbx6+)", "Somites / Dermomyotome" ,"Somites / Sclerotome"  ))

test_pair <- c("HAP")
subset <- subset(Asmb_somites, condition == test_pair[1])
subset <- subset(subset, predicted.celltype_updated %in% c("Mesodermal progenitors (Tbx6+)", "Somites / Dermomyotome" ,"Somites / Sclerotome"  ))
subset <- ScaleData(subset, features= markers, verbose = FALSE)
Idents(subset) <- "predicted.celltype_updated"

scale_mat <- GetAssayData(subset, assay = "RNA", layer = "scale.data")
present <- intersect(markers, rownames(scale_mat))
missing <- setdiff(markers, rownames(scale_mat))

# zero matrix for missing genes
if (length(missing) > 0) {
  zero_mat <- matrix(
    0,
    nrow = length(missing),
    ncol = ncol(scale_mat),
    dimnames = list(missing, colnames(scale_mat))
  )
  subset_data <- rbind(scale_mat[present, , drop = FALSE], zero_mat)
} else {
  subset_data <- scale_mat[present, , drop = FALSE]
}

gene_order <- unique(unlist(marker_groups, use.names = FALSE))
subset_data <- subset_data[gene_order[gene_order %in% rownames(subset_data)], , drop = FALSE]
subset_meta <- subset@meta.data
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell") 

df <- as.data.frame(t(as.matrix(subset_data)))
df$cell <- rownames(df)
subset_meta <- subset@meta.data %>%
  rownames_to_column(var = "cell")  
df <- df %>%
  pivot_longer(cols = -cell, names_to = "gene", values_to = "expression") %>%
  left_join(subset_meta, by = "cell") 

# mean expression per group
df_summary <- df %>%
  group_by(gene,  condition, predicted.celltype_updated) %>%
  summarize(
    avg_expression = mean(expression), 
    num_expressing = sum(expression > 0),  
    total_cells = n(),  
    pct_expressed = (num_expressing / total_cells) * 100,  
    .groups = "drop"
  )

subset <- df_summary
gene_order <- unlist(marker_groups, use.names = FALSE) 
subset$gene <- factor(subset$gene, levels = gene_order)

# DotPlot
filename <- paste0("images/scRNAseq/OMG/marker_expression_Somites_",test_pair[1], ".pdf")
pdf(file = filename,height=5, width=10)
ggplot(subset, aes(x = gene, y = predicted.celltype_updated, size = pct_expressed)) +
  geom_point(shape = 21, color = "black", stroke = 0.5,  aes(fill = avg_expression)) +
  scale_fill_gradient2(low = "#2166ac", mid = "gray93", high = "#b2182b", midpoint = 0) +
  scale_size(range = c(1, 8)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    text = element_text(color = "black"),              
    axis.text = element_text(color = "black"),         
    axis.title = element_text(color = "black"),        
    strip.text = element_text(color = "black"),        
    legend.text = element_text(color = "black"),       
    legend.title = element_text(color = "black")       
  ) +
  labs(x = "", y = "", fill = "Avg. Expression", size = "Percent Expressed") + coord_fixed(ratio = 9/7) 
dev.off()
