
# Cross-species staging against human embryo atlases (Xu 2023 + Zeng 2023)
#
# Input : data/scRNAseq/human/GSE157329_* (Xu et al. 2023 raw matrix + gene/
#         cell annotation), data/scRNAseq/human/GSE155121_Human.rds (Zeng et al.
#         2023), and data/scRNAseq/Asmb_OMG_filtered.rds (from 04)
#
# Steps :
#   Xu 2023 
#     1. Build a Seurat object, QC, normalize/PCA/UMAP
#     2. Human->mouse 1:1 orthologs (homologene) and subset to genes shared with the gastruloid data
#     3. RPCA-integrate the v2/v3 kit batches
#     4. Transfer Carnegie-stage labels to the gastruloids (equal-per-stage
#        subsample), per-condition heatmaps, stage dot plots, and 
#        score-weighted expected stage per cell type
#   Zeng 2023 
#     5. QC + per-sample DoubletFinder (pK swept per week_stage, homotypic
#        correction, size-scaled expected rate capped at 6%)
#     6. Human->mouse orthologs via explicit matrix extraction and CreateAssay5Object
#     7. Transfer post-conceptual-week labels (with / without RPCA integration,
#        with / without per-sample suffix); heatmaps, dot plots, weighted stage
#
# Output: data/scRNAseq/human/{Human_Zeng_sub, Asmb_Zeng_sub,
#         Human_Zeng_eqStage_sub}.rds; figures in images/scRNAseq/Human/
#         and images/scRNAseq/Human/Zeng_2023/

###############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(DoubletFinder)
library(reticulate)
library(ggplot2)
library(homologene)
library(tidyr)
set.seed(42)
source("scripts/OMG_colors.r")

# rename predicted.id / predicted.id.score to avoid column conflicts when adding all three
rename_transfer <- function(df, tag) {
  colnames(df)[colnames(df) == "predicted.id"]       <- paste0("pred.", tag)
  colnames(df)[colnames(df) == "predicted.id.score"] <- paste0("pred.", tag, ".score")
  df
}

# for each condition plot heatmap showing predicted human labels per cell type in Asmb
make_heatmaps <- function(seurat_obj, pred_col, y_label, out_pdf) {
  df <- seurat_obj@meta.data[, c("condition", "predicted.celltype_updated", pred_col)]
  colnames(df)[3] <- "pred_label"

  prop_df <- df %>%
    group_by(condition, predicted.celltype_updated, pred_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(condition, predicted.celltype_updated) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

  conditions <- unique(prop_df$condition)

  plot_list <- lapply(conditions, function(cond) {
    ggplot(
      prop_df[prop_df$condition == cond, ],
      aes(x = predicted.celltype_updated, y = pred_label, fill = proportion)
    ) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient(low = "white", high = "#2166AC", name = "Proportion") +
      labs(
        title = cond,
        x = "Predicted cell type (Asmb)",
        y = y_label
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold")
      )
  })

  pdf(out_pdf, width = 14, height = 8)
  for (p in plot_list) print(p)
  dev.off()
}

data_dir <- "data/scRNAseq/human"

# read the three files
counts   <- readMM(gzcon(file(file.path(data_dir, "GSE157329_raw_counts.mtx.gz"), "rb")))
cell_ann <- read.table(file.path(data_dir, "GSE157329_cell_annotate.txt.gz"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_ann <- read.table(file.path(data_dir, "GSE157329_gene_annotate.txt.gz"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# MTX 
rownames(counts) <- make.unique(gene_ann$gene_short_name)
colnames(counts) <- cell_ann$cell_id

# create Seurat object and attach cell metadata
human <- CreateSeuratObject(counts  = counts,
                            project = "GSE157329",
                            min.cells = 3, min.features = 200)

# add all metadata columns from cell_annotate
meta <- cell_ann %>% column_to_rownames("cell_id")
human <- AddMetaData(human, metadata = meta)

# save the Seurat object for downstream analysis
#saveRDS(human, file.path(data_dir, "Human.rds"))
#human <- readRDS(file.path(data_dir, "Human.rds"))
# mito content
human[["percent.mt"]] <- PercentageFeatureSet(human, pattern = "^MT-")
Asmb <- readRDS("data/scRNAseq/Asmb_OMG_filtered.rds")

human <- NormalizeData(human)
human <- FindVariableFeatures(human, selection.method = "vst", nfeatures = 2000)
human <- ScaleData(human)
human <- RunPCA(human, features = VariableFeatures(object = human))
human <- RunUMAP(human, dims = 1:30, return.model = TRUE)

# translate human gene symbols to mouse orthologs via homologene 
orthologs <- human2mouse(rownames(human), db = homologeneData)
colnames(orthologs) <- c("human_gene", "mouse_gene")

# keep only 1-to-1 mappings 
keep <- !duplicated(orthologs$human_gene) &
  !duplicated(orthologs$mouse_gene)
orthologs <- orthologs[keep, ]

# subset and rename the human object to mouse gene symbols
human_genes_keep <- intersect(rownames(human), orthologs$human_gene)
orthologs <- orthologs[orthologs$human_gene %in% human_genes_keep, ]
human_sub <- human[human_genes_keep, ]
mouse_names <- orthologs$mouse_gene[
  match(human_genes_keep, orthologs$human_gene)
]
rownames(human_sub[["RNA"]]) <- mouse_names

# find shared genes with Asmb and subset both objects
shared_genes <- intersect(rownames(human_sub), rownames(Asmb))
human_sub <- human_sub[shared_genes, ]
asmb_sub  <- Asmb[shared_genes, ]

# create batch variable matching paper (v2: Emb01/02/05/06/07, v3: Emb03/04)
human_sub$kit_batch <- ifelse(
  human_sub$embryo %in% c("emb3", "emb4"), "v3", "v2"
)

# split RNA layers by batch
human_sub[["RNA"]] <- split(human_sub[["RNA"]], f = human_sub$kit_batch)

human_sub <- FindVariableFeatures(human_sub, selection.method = "vst", nfeatures = 2000)
human_sub <- ScaleData(human_sub)
human_sub <- RunPCA(human_sub, features = VariableFeatures(human_sub))

human_sub <- IntegrateLayers(
  human_sub,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  dims           = 1:30
)
human_sub <- JoinLayers(human_sub)
human_sub <- RunUMAP(human_sub, reduction = "integrated.rpca",
                     dims = 1:30, return.model = TRUE)

out_dir <- "images/scRNAseq/Human"
# plot human umap 
pdf(file.path(out_dir, "Human_UMAP.pdf"), width = 8, height = 6)
DimPlot(human_sub, reduction = "umap", group.by = "stage") +
  labs(title = "Human scRNA-seq UMAP colored by stage") 
dev.off()

pdf(file.path(out_dir, "Human_UMAP_by_dev_system.pdf"), width = 10, height = 8)
DimPlot(human_sub, reduction = "umap", group.by = "developmental.system") +
  labs(title = "Human scRNA-seq UMAP colored by developmental system")
dev.off()

pdf(file.path(out_dir, "Human_UMAP_by_embryo.pdf"), width = 10, height = 8)
DimPlot(human_sub, reduction = "umap", group.by = "embryo") +
  labs(title = "Human scRNA-seq UMAP colored by embryo")
dev.off()


# subsample human to same amount of cells per stage 
human_eq_stage_cells <- human_sub@meta.data %>%
  # cell_id = colnames(human_sub)
  mutate(cell_id = rownames(human_sub@meta.data)) %>%
  subset(stage %in% c("CS12", "CS13-14", "CS15-16")) %>%
  group_by(stage) %>%
  sample_n(25250) %>%
  ungroup()

human_eq_stage <- subset(human_sub , cells=human_eq_stage_cells$cell_id)
human_eq_stage <- NormalizeData(human_eq_stage)
human_eq_stage <- FindVariableFeatures(human_eq_stage, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(human_eq_stage)
human_eq_stage <- ScaleData(human_eq_stage, features = all.genes)
human_eq_stage <- RunPCA(human_eq_stage)
# project Asmb onto human
anchors <- FindTransferAnchors(
  reference           = human_eq_stage,
  query               = asmb_sub,
  reference.reduction = "integrated.rpca",
  dims                = 1:30
)
pred_stage  <- TransferData(anchorset = anchors, refdata = human_eq_stage$stage,dims = 1:30)
pred_stage  <- rename_transfer(pred_stage,  "stage_eq")
Asmb <- AddMetaData(Asmb, pred_stage)

make_heatmaps(Asmb, "pred.stage_eq", "Stage (human_stage_eq)",
              file.path(out_dir, "Asmb_stage_per_condition_eq.pdf"))

table(human_eq_stage$stage)
table(Asmb$pred.stage_eq)
df <- Asmb@meta.data


# DotPlot to show cell fraction per predicted stage 
dotplot_data <- df %>%
  dplyr::group_by(condition, pred.stage_eq) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%  # normalization within each condition
  dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) 
print(dotplot_data)

unique(human$stage)
all_stages <- c("CS12" , "CS13-14" , "CS15-16")
dotplot_data$pred.stage <- factor(dotplot_data$pred.stage_eq, levels = all_stages, ordered = TRUE)

dotplot_data <- dotplot_data %>%
    tidyr::complete(pred.stage = all_stages,
                    fill = list(condition = unique(dotplot_data$condition), cell_count = NA, percentage = NA, pred.stage = NA))
dotplot_data$pred.stage <- factor(dotplot_data$pred.stage, levels = all_stages, ordered = TRUE)
  
plot <- ggplot(dotplot_data, aes(x = pred.stage, y = condition, size = percentage)) +
    geom_point(alpha = 0.8, na.rm = TRUE) + 
    scale_size_continuous(range = c(2, 10)) +  # No dots for 0% values
    scale_x_discrete(expand = c(0.1, 0.1)) +
     coord_cartesian(clip = "off") +
    labs(x = "Predicted human stage (CS)", y = "", size = "Cell fraction (%)", title = "") +
    theme_minimal() +
    theme_classic() +
    theme(
      # rotate x-axis labels for better readability
      axis.text.x = element_text(size = 14, color = "black", angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),  
      axis.title = element_text(size = 16),  
      legend.text = element_text(size = 12),  
      legend.title = element_text(size = 14), aspect.ratio = 0.8
)
filename <- paste0("images/scRNAseq/Human/Dotplot_stage.pdf")
ggsave(filename, width = 6, height = 3, plot)


# do boxplot for all conditions
for(cond in unique(df$condition)){
  cells_cond <- rownames(subset(df, condition == cond))
  pred_week_scores_cond <- pred_stage[cells_cond,]
  
  # extract only per-stage score columns (exclude character and max columns)
  score_cols <- grep("^prediction\\.score\\.CS", colnames(pred_week_scores_cond), value = TRUE)
  score_mat  <- as.matrix(pred_week_scores_cond[, score_cols])
  
  # CS number from column name 
  cs_map <- c("prediction.score.CS12" = 12,
            "prediction.score.CS13.14" = 13.5,
            "prediction.score.CS15.16" = 15.5)
  w <- cs_map[colnames(score_mat)]
  print(head(score_mat))
  # Weighted sum (expected week) per cell
  expected_day <- (score_mat %*% w) / rowSums(score_mat)
  
  df_w <- data.frame(
    cell         = rownames(pred_week_scores_cond),
    expected_day = as.numeric(expected_day),
    row.names    = rownames(pred_week_scores_cond)
  )
  
  # add information about cell type
  df_w$cell_type <- Asmb$predicted.celltype_updated[rownames(pred_week_scores_cond)]
  
  df_w$cell_type <- factor(df_w$cell_type,
                           levels = names(cell_type_colored))
  
  pdf(paste0("images/scRNAseq/Human/Boxplot_stage_weighted_sum_", cond, "_per_celltype.pdf"),
      width = 14, height = 6)
  p <- ggplot(df_w, aes(x = cell_type, y = expected_day,  fill = cell_type)) +
    geom_boxplot(outlier.alpha = 0.3) +
    labs(y = "Predicted human stage (PCW)", x = "", title = paste0("Condition: ", cond)) +
     scale_fill_manual(values = cell_type_colored) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      plot.title = element_text(size=18)
    )
    print(p)
  dev.off()
}


#############################################################################################
#############################################################################################
# Zeng et al. 2023 
human_Zeng <- readRDS("data/scRNAseq/human/GSE155121_Human.rds")
human <- human_Zeng
Asmb <- readRDS("data/scRNAseq/Asmb_OMG_filtered.rds")

# QC human
# get nFeature_RNA and nCount_RNA
human$nFeature_RNA <- colSums(GetAssayData(human, assay = "RNA", layer = "counts") > 0)
human$nCount_RNA   <- colSums(GetAssayData(human, assay = "RNA", layer = "counts"))
human$percent.mt <- PercentageFeatureSet(human, pattern = "^MT-")
# without points
out_dir <- "images/scRNAseq/Human/Zeng_2023"
pdf(file.path(out_dir, "Human_Zeng_QC.pdf"), width = 12, height = 4)
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
dev.off()

# as in Zeng et al. 2023
human <- subset(human, subset = nCount_RNA > 200 & nCount_RNA< 50000  & nFeature_RNA>800 & nFeature_RNA < 6000  & percent.mt < 20)

human[["RNA"]]@meta.features <- data.frame(row.names = rownames(human))
human <- NormalizeData(human)
human <- FindVariableFeatures(human)
human <- ScaleData(human)
human <- RunPCA(human)
human <- RunUMAP(human, dims = 1:30)
human <- FindNeighbors(human, dims = 1:30)
human <- FindClusters(human)

# Sweep pK (optimal neighbourhood size)
sweep_res <- paramSweep(human, PCs = 1:30, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
pdf(file.path(out_dir, "DoubletFinder_pK_sweep.pdf"))
bcmvn <- find.pK(sweep_stats)
dev.off()
optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Estimate homotypic doublet proportion
annotations <- human$seurat_clusters
homotypic_prop <- modelHomotypic(annotations)
nExp <- round(0.06 * ncol(human) * (1 - homotypic_prop))  # 6% expected doublet rate

human_v4 <- human
human_v4[["RNA"]] <- as(human_v4[["RNA"]], "Assay")

results <- list()
for (sample in unique(human_v4$week_stage)) {
  cat("Processing sample:", sample, "\n")
  
  sub <- subset(human_v4, week_stage == sample)
  
  # Re-run preprocessing on subset
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:30, verbose = FALSE)
  sub <- FindClusters(sub, verbose = FALSE)
  
  # pK sweep per sample
  sweep_res <- paramSweep(sub, PCs = 1:30, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  pdf(tempfile(fileext = ".pdf"))
  bcmvn_sub <- find.pK(sweep_stats)
  dev.off()
  pk_sub <- as.numeric(as.character(
    bcmvn_sub$pK[which.max(bcmvn_sub$BCmetric)]
  ))
  
  # Doublet rate based on per-sample cell count
  n_cells <- ncol(sub)
  rate <- min(0.06, 0.008 * (n_cells / 1000))
  annotations_sub <- sub$seurat_clusters
  homotypic_sub <- modelHomotypic(annotations_sub)
  nExp_sub <- round(rate * n_cells * (1 - homotypic_sub))
  
  sub <- doubletFinder(sub, PCs = 1:30, pN = 0.25, pK = pk_sub,
                       nExp = nExp_sub, reuse.pANN = NULL, sct = FALSE)
  
  df_col <- grep("DF.classifications", colnames(sub@meta.data), value = TRUE)
  results[[sample]] <- sub@meta.data[, df_col, drop = FALSE]
}

# Combine results back into main object
all_classifications <- do.call(rbind, lapply(names(results), function(s) {
  df <- results[[s]]
  colnames(df) <- "doublet_class"
  df
}))

human_v4$doublet_class <- all_classifications[rownames(human_v4@meta.data), "doublet_class"]
human$doublet_class <- human_v4$doublet_class

# Check how many doublets
table(human$doublet_class)

# Filter
human <- subset(human, doublet_class == "Singlet")
cat("Cells after doublet removal:", ncol(human), "\n")
# Cells after doublet removal: 417126 
human$all_cells <- "all"

out_dir <- "images/scRNAseq/Human/Zeng_2023"
pdf(file.path(out_dir, "Human_Zeng_after_QC.pdf"), width = 10, height = 4)
VlnPlot(human, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="all_cells",ncol = 3, pt.size = 0) 
dev.off()

dim(human)
head(human[["RNA"]])

# translate human gene symbols to mouse orthologs via homologene (local DB)
orthologs <- human2mouse(rownames(human), db = homologeneData)
colnames(orthologs) <- c("human_gene", "mouse_gene")

# keep only 1-to-1 mappings 
keep <- !duplicated(orthologs$human_gene) &
  !duplicated(orthologs$mouse_gene)
orthologs <- orthologs[keep, ]

# subset and rename the human object to mouse gene symbols
human_genes_keep <- intersect(rownames(human), orthologs$human_gene)
orthologs <- orthologs[orthologs$human_gene %in% human_genes_keep, ]
human_sub <- human[human_genes_keep, ]
mouse_names <- orthologs$mouse_gene[
  match(human_genes_keep, orthologs$human_gene)
]

counts_mat <- GetAssayData(human_sub, assay = "RNA", layer = "counts")
data_mat   <- GetAssayData(human_sub, assay = "RNA", layer = "data")

# Rename rows to mouse symbols
rownames(counts_mat) <- mouse_names
rownames(data_mat)   <- mouse_names

# Build v5 assay with mouse names
new_assay <- CreateAssay5Object(counts = counts_mat)
new_assay <- SetAssayData(new_assay, layer = "data", new.data = data_mat)
human_sub[["RNA"]] <- new_assay

head(rownames(human_sub[["RNA"]]))
length(rownames(human_sub[["RNA"]]))

# find shared genes with Asmb and subset both objects
shared_genes <- intersect(rownames(human_sub), rownames(Asmb))
human_sub <- human_sub[shared_genes, ]
asmb_sub  <- Asmb[shared_genes, ]

# # integrate week_stage batches in human_sub with rpca
# # split RNA layers by batch so IntegrateLayers can operate per-batch
# human_sub[["RNA"]] <- split(human_sub[["RNA"]], f = human_sub$week_stage)

human_sub <- NormalizeData(human_sub)
human_sub <- FindVariableFeatures(human_sub, selection.method = "vst", nfeatures = 2000)
human_sub <- ScaleData(human_sub)
human_sub <- RunPCA(human_sub, features = VariableFeatures(human_sub))

# save human_sub for downstream analysis
saveRDS(human_sub, file.path("data/scRNAseq/human/Human_Zeng_sub.rds"))

out_dir <- "images/scRNAseq/Human/Zeng_2023"

asmb_sub <- NormalizeData(asmb_sub)
asmb_sub <- FindVariableFeatures(asmb_sub, selection.method = "vst", nfeatures = 2000)
asmb_sub <- ScaleData(asmb_sub)
asmb_sub <- RunPCA(asmb_sub, features = VariableFeatures(asmb_sub))

table(human_sub$week_stage)

human_eq_week_cells <- human_sub@meta.data %>%
  mutate(cell_id = rownames(human_sub@meta.data),
         week    = sub("-.*$", "", week_stage)) %>%   # W4-1/W4-2/W4-3 -> W4
  group_by(week) %>%
  sample_n(11500) %>%                                  # equal cells PER WEEK
  ungroup()

human_sub$week_stage_woSample <- gsub("-.*", "", human_sub$week_stage)
human_eq_stage <- subset(human_sub , cells=human_eq_week_cells$cell_id)
table(human_eq_stage$week_stage)

human_eq_stage <- NormalizeData(human_eq_stage)
human_eq_stage <- FindVariableFeatures(human_eq_stage, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(human_eq_stage)
human_eq_stage <- ScaleData(human_eq_stage)
human_eq_stage <- RunPCA(human_eq_stage)

anchors <- FindTransferAnchors(
  reference           = human_eq_stage,
  query               = asmb_sub,
  reference.reduction = "pca",
  dims                = 1:50
)

pred_week  <- TransferData(anchorset = anchors, refdata = human_eq_stage$week_stage,  dims = 1:50)
pred_week_woSample  <- TransferData(anchorset = anchors, refdata = human_eq_stage$week_stage_woSample,  dims = 1:50)

unique(pred_week$predicted.id)

pred_week  <- rename_transfer(pred_week,  "week_stage")
pred_week_woSample  <- rename_transfer(pred_week_woSample,  "week_stage_woSample")
Asmb <- AddMetaData(Asmb, pred_week)
Asmb <- AddMetaData(Asmb, pred_week_woSample)

make_heatmaps(Asmb, "pred.week_stage", "Week stage (human_week_stage)",
              file.path(out_dir, "Asmb_week_stage_per_condition_woIntegration.pdf"))

make_heatmaps(Asmb, "pred.week_stage_woSample", "Week stage (human_week_stage_woSample)",
              file.path(out_dir, "Asmb_week_stage_per_condition_woSample_woIntegration.pdf"))

df <- Asmb@meta.data


## Staging, dotplot for HAPs
df <- Asmb@meta.data
haps <- subset(df, condition == "HAP")

# DotPlot to show cell fraction per predicted stage 
dotplot_data <- haps %>%
  dplyr::group_by(condition, pred.week_stage_woSample) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%  # normalization within each condition
  dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) 
print(dotplot_data)


unique(human_sub$week_stage_woSample)
all_stages <- c("W3" , "W4" , "W5" , "W6" , "W7" , "W8" , "W9" , "W12")
# remove -any number from labels in pred.week_stage
dotplot_data$predicted.week <- gsub("-\\d+", "", dotplot_data$pred.week_stage_woSample)
# ensure order of stages in plot  
dotplot_data$predicted.week <- factor(dotplot_data$predicted.week, levels = all_stages, ordered = TRUE)


  dotplot_data <- dotplot_data %>%
    tidyr::complete(predicted.week = all_stages,
                    fill = list(condition = "HAP", cell_count = NA, percentage = NA, pred.week_stage = NA))
  
  dotplot_data$predicted.week <- factor(dotplot_data$predicted.week, levels = all_stages, ordered = TRUE)
  

  plot <- ggplot(dotplot_data, aes(x = predicted.week, y = condition, size = percentage)) +
    geom_point(alpha = 0.8, na.rm = TRUE) + 
    scale_size_continuous(range = c(2, 10)) +  # No dots for 0% values
    scale_x_discrete(expand = c(0.1, 0.1)) +
     coord_cartesian(clip = "off") +
    labs(x = "Predicted human stage (PCW)", y = "", size = "Cell fraction (%)", title = "") +
    theme_minimal() +
    theme_classic() +
    theme(
      # rotate x-axis labels for better readability
      axis.text.x = element_text(size = 14, color = "black", angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),  
      axis.title = element_text(size = 16),  
      legend.text = element_text(size = 12),  
      legend.title = element_text(size = 14), aspect.ratio = 0.5
    )
  filename <- paste0("images/scRNAseq/Human/Zeng_2023/Dotplot_stage_HAP.pdf")
  ggsave(filename, width = 7, height = 2.8, plot)

# for all conditions
# DotPlot to show cell fraction per predicted stage 
dotplot_data <- df %>%
  dplyr::group_by(condition, pred.week_stage_woSample) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%  # normalization within each condition
  dplyr::mutate(percentage = (cell_count / sum(cell_count)) * 100) 
print(dotplot_data)


unique(human_sub$week_stage_woSample)
all_stages <- c("W3" , "W4" , "W5" , "W6" , "W7" , "W8" , "W9" , "W12")
# remove -any number from labels in pred.week_stage
dotplot_data$predicted.week <- gsub("-\\d+", "", dotplot_data$pred.week_stage_woSample)
# ensure order of stages in plot  
dotplot_data$predicted.week <- factor(dotplot_data$predicted.week, levels = all_stages, ordered = TRUE)


  dotplot_data <- dotplot_data %>%
    tidyr::complete(predicted.week = all_stages,
                    fill = list(condition = unique(dotplot_data$condition), cell_count = NA, percentage = NA, pred.week_stage = NA))
  
  dotplot_data$predicted.week <- factor(dotplot_data$predicted.week, levels = all_stages, ordered = TRUE)
  

  plot <- ggplot(dotplot_data, aes(x = predicted.week, y = condition, size = percentage)) +
    geom_point(alpha = 0.8, na.rm = TRUE) + 
    scale_size_continuous(range = c(2, 10)) +  # No dots for 0% values
    scale_x_discrete(expand = c(0.1, 0.1)) +
     coord_cartesian(clip = "off") +
    labs(x = "Predicted human stage (PCW)", y = "", size = "Cell fraction (%)", title = "") +
    theme_minimal() +
    theme_classic() +
    theme(
      # rotate x-axis labels for better readability
      axis.text.x = element_text(size = 14, color = "black", angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),  
      axis.title = element_text(size = 16),  
      legend.text = element_text(size = 12),  
      legend.title = element_text(size = 14), aspect.ratio = 0.5
    )
  filename <- paste0("images/scRNAseq/Human/Zeng_2023/Dotplot_stage.pdf")
  ggsave(filename, width = 7, height = 2.8, plot)


# Staging per cell type 

hapa_cells <- rownames(haps)
pred_week_scores <- pred_week_woSample[hapa_cells,]


# extract only per-stage score columns (exclude character and max columns)
score_cols <- grep("^prediction\\.score\\.W", colnames(pred_week_scores), value = TRUE)
score_mat  <- as.matrix(pred_week_scores[, score_cols])

# weights: week number from column name (e.g. prediction.score.W4 -> 4)
w <- as.numeric(sub("prediction\\.score\\.W", "", score_cols))

# 3) Weighted sum (expected week) per cell
expected_day <- (score_mat %*% w) / rowSums(score_mat)

df_w <- data.frame(
  cell         = rownames(pred_week_scores),
  expected_day = as.numeric(expected_day),
  row.names    = rownames(pred_week_scores)
)

# add information about cell type
df_w$cell_type <- haps[rownames(pred_week_scores), "predicted.celltype_updated"]

df_w$cell_type <- factor(df_w$cell_type,
                         levels = names(cell_type_colored))

pdf("images/scRNAseq/Human/Zeng_2023/Boxplot_stage_weighted_sum_HAPs_per_celltype.pdf",
    width = 12, height = 6)
ggplot(df_w, aes(x = cell_type, y = expected_day,  fill = cell_type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(y = "Predicted human stage (PCW)", x = "") +
   scale_fill_manual(values = cell_type_colored) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()



# do boxplot for all conditions
for(cond in unique(df$condition)){
  cells_cond <- rownames(subset(df, condition == cond))
  pred_week_scores_cond <- pred_week_woSample[cells_cond,]
  
  # extract only per-stage score columns (exclude character and max columns)
  score_cols <- grep("^prediction\\.score\\.W", colnames(pred_week_scores_cond), value = TRUE)
  score_mat  <- as.matrix(pred_week_scores_cond[, score_cols])
  
  # weights: week number from column name (e.g. prediction.score.W4 -> 4)
  w <- as.numeric(sub("prediction\\.score\\.W", "", score_cols))
  
  # Weighted sum (expected week) per cell
  expected_day <- (score_mat %*% w) / rowSums(score_mat)
  
  df_w <- data.frame(
    cell         = rownames(pred_week_scores_cond),
    expected_day = as.numeric(expected_day),
    row.names    = rownames(pred_week_scores_cond)
  )
  
  # add information about cell type
  df_w$cell_type <- asmb_sub$predicted.celltype_updated[rownames(pred_week_scores_cond)]
  
  df_w$cell_type <- factor(df_w$cell_type,
                           levels = names(cell_type_colored))
  
  pdf(paste0("images/scRNAseq/Human/Zeng_2023/Boxplot_stage_weighted_sum_", cond, "_per_celltype.pdf"),
      width = 12, height = 6)
  p <- ggplot(df_w, aes(x = cell_type, y = expected_day,  fill = cell_type)) +
    geom_boxplot(outlier.alpha = 0.3) +
    labs(y = "Predicted human stage (PCW)", x = "", title = paste0("Condition: ", cond)) +
     scale_fill_manual(values = cell_type_colored) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      plot.title = element_text(size=18)
    )
    print(p)
  dev.off()
}

# save human_sub
saveRDS(human_sub, file.path("data/scRNAseq/human/Human_Zeng_sub.rds"))
saveRDS(asmb_sub, file.path("data/scRNAseq/human/Asmb_Zeng_sub.rds"))
saveRDS(human_eq_stage, file.path("data/scRNAseq/human/Human_Zeng_eqStage_sub.rds"))
