#
# Build the OMG ("Ontogeny of Mouse, Graphed") reference (E8.0-E10.0)
#
# Input : Public OMG/JAX scRNA-seq downloads (see URLs below),
#         cell metadata (df_cell.rds), gene annotation, and count matrices
#         for run_4, run_15, run_17_sub1, run_17_sub2. Expected in data/OMG/.
#
# Steps : 1. Load cell metadata; restrict to embryonic days E8.0-E10.0
#         2. Load per-run count matrices, attach gene/cell IDs, subset to the
#            selected cells, and build per-run Seurat objects
#         3. QC filter (nFeature 1500-3000, nCount 2500-6000)
#         4. Downsample each stage to 15,000 cells and merge into one object
#         5. Map Ensembl IDs to gene symbols, dropping non-unique symbols, and
#            rebuild the assay for consistency
#         6. Normalize, scale, PCA, UMAP; QC and overview plots
#
# Output: data/OMG/OMG_E8_E10_merge.rds  (input to 02_annotate_OMG_reference.R)
#         reference overview figures in images/OMG/
#
# Note  : This is the unrestricted E8.0-E10.0 reference. Restriction to
#         E8.0-E9.75 and annotation refinement happen in 02.

###############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

library(Seurat)
library(future)
library(future.apply)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrastr)
source("data/OMG/JAX_color_code.R")

##############################################################################
# Load cell metadata
# download link: https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/meta_data/df_cell.rds
pd_all <- readRDS("data/OMG/df_cell.rds")

# Select specific days (E8.0 to E10.0)
pd <- subset(pd_all, day %in% c('E8.0-E8.5','E8.75','E9.0','E9.25','E9.5','E9.75','E10.0'))

# Download individual datasets (run_4 to run_28)
# Example URLs:
# https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/mtx/gene_count.run_4.mtx.gz
# https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/mtx/cell_annotation.run_4.csv.gz
# https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/mtx/gene_annotation.csv.gz

unique(gsub('_P.*','',pd$cell_id))
# [1] "run_4"  "run_15" "run_17"

# Load data and make Seurat object
df_gene <- read.csv("data/OMG/gene_annotation.csv.gz")

## run4
gene_count <- Matrix::readMM("data/OMG/gene_count.run_4.mtx.gz")
df_cell <- read.csv("data/OMG/cell_annotation.run_4.csv.gz")
### Assign row and column names
rownames(gene_count) <- df_gene$gene_ID
colnames(gene_count) <- df_cell$cell_id
### Subset gene count matrix
gene_count <- gene_count[, colnames(gene_count) %in% pd$cell_id, drop = FALSE]
pd.v <- pd
rownames(pd.v) = as.vector(pd.v$cell_id)
pd.v = pd.v[colnames(gene_count),]
WT.run4 <- CreateSeuratObject(gene_count, meta.data = pd.v)

## run15
gene_count <- Matrix::readMM("data/OMG/gene_count.run_15.mtx.gz")
df_cell <- read.csv("data/OMG/cell_annotation.run_15.csv.gz")
### Assign row and column names
rownames(gene_count) <- df_gene$gene_ID
colnames(gene_count) <- df_cell$cell_id
### Subset gene count matrix
gene_count <- gene_count[, colnames(gene_count) %in% pd$cell_id, drop = FALSE]
pd.v <- pd
rownames(pd.v) = as.vector(pd.v$cell_id)
pd.v = pd.v[colnames(gene_count),]
WT.run15 <- CreateSeuratObject(gene_count, meta.data = pd.v)

## run17_sub1
gene_count <- Matrix::readMM("data/OMG/gene_count.run_17_sub1.mtx.gz")
df_cell <- read.csv("data/OMG/cell_annotation.run_17_sub1.csv.gz")
### Assign row and column names
rownames(gene_count) <- df_gene$gene_ID
colnames(gene_count) <- df_cell$cell_id
### Subset gene count matrix
gene_count <- gene_count[, colnames(gene_count) %in% pd$cell_id, drop = FALSE]
pd.v <- pd
rownames(pd.v) = as.vector(pd.v$cell_id)
pd.v = pd.v[colnames(gene_count),]
WT.run17.1 <- CreateSeuratObject(gene_count, meta.data = pd.v)

## run17_sub2
gene_count <- Matrix::readMM("data/OMG/gene_count.run_17_sub2.mtx.gz")
df_cell = read.csv("data/OMG/cell_annotation.run_17_sub2.csv.gz")
### Assign row and column names
rownames(gene_count) <- df_gene$gene_ID
colnames(gene_count) <- df_cell$cell_id
### Subset gene count matrix
gene_count <- gene_count[, colnames(gene_count) %in% pd$cell_id, drop = FALSE]
pd.v <- pd
rownames(pd.v) = as.vector(pd.v$cell_id)
pd.v = pd.v[colnames(gene_count),]
WT.run17.2 <- CreateSeuratObject(gene_count, meta.data = pd.v)


# QC and filtering
WT.run4.x <- subset(WT.run4, subset = nFeature_RNA > 1500 & nFeature_RNA < 3000 & nCount_RNA > 2500 & nCount_RNA < 6000)
WT.run15.x <- subset(WT.run15, subset = nFeature_RNA > 1500 & nFeature_RNA < 3000 & nCount_RNA > 2500 & nCount_RNA < 6000)
WT.run17.1.x <- subset(WT.run17.1, subset = nFeature_RNA > 1500 & nFeature_RNA < 3000 & nCount_RNA > 2500 & nCount_RNA < 6000)
WT.run17.2.x <- subset(WT.run17.2, subset = nFeature_RNA > 1500 & nFeature_RNA < 3000 & nCount_RNA > 2500 & nCount_RNA < 6000)

# Downsample for equal day distribution and merge Seurat objects
WT.E8.5 <- subset(WT.run4.x, downsample = 15000)
WT.E8.75 <- subset(subset(WT.run15.x, day=='E8.75'), downsample = 15000)
WT.E9.0 <- subset(subset(WT.run15.x, day=='E9.0'), downsample = 15000)
WT.E9.25 <- subset(subset(WT.run15.x, day=='E9.25'), downsample = 15000)
WT.E9.5 <- subset(subset(WT.run15.x, day=='E9.5'), downsample = 15000)
WT.E9.75 <- subset(subset(merge(x=WT.run17.1.x, y=WT.run17.2.x), day=='E9.75'), downsample = 15000)
WT.E10.0 <- subset(subset(merge(x=WT.run17.1.x, y=WT.run17.2.x), day=='E10.0'), downsample = 15000)

OMG <- merge(x=WT.E8.5, y=c(WT.E8.75, WT.E9.0, WT.E9.25, WT.E9.5, WT.E9.75, WT.E10.0))
OMG@meta.data$day <- factor(OMG@meta.data$day, levels = c("E8.0-E8.5","E8.75", "E9.0","E9.25","E9.5","E9.75", "E10.0"))

pdf('images/OMG/all_states_cell_type_distribution.pdf', height=120, width=21)
ggplot(OMG@meta.data, aes(x=day, fill=celltype_update)) + geom_bar() + theme_classic() + facet_grid(celltype_update~., scales='free_y') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

quantile(OMG@meta.data$nFeature_RNA)
#   0%  25%  50%  75% 100%
# 1501 1928 2216 2546 2999
quantile(OMG@meta.data$nCount_RNA)
#   0%  25%  50%  75% 100%
# 2501 3070 3731 4570 5999
OMG
# An object of class Seurat
# 49585 features across 105000 samples within 1 assay
# Active assay: RNA (49585 features, 0 variable features)
#  8 layers present: counts.1, counts.2, counts.3, counts.4, counts.5, counts.1.6, counts.2.6, counts.2.7

OMG <- JoinLayers(OMG)
OMG_Ensemble <- OMG

# Rename genes
### Get the counts and strip Ensembl version suffixes
assay_name <- DefaultAssay(OMG)
counts <- GetAssayData(OMG, assay = assay_name, slot = "counts")
ens_ids <- gsub("\\..*$", "", rownames(counts))

### Map ID -> name
symbols <- df_gene[,c('gene_ID','gene_short_name')]
rownames(symbols) <- df_gene$gene_ID
symbols <- symbols[ens_ids,]
symbols_unique <- subset(symbols, !gene_short_name %in% symbols$gene_short_name[duplicated(symbols$gene_short_name)])

### Remove genes that are non-unique
counts_unique <- counts[which(rownames(counts) %in% symbols_unique$gene_ID),]
rownames(counts_unique) <- symbols_unique[rownames(counts_unique),]$gene_short_name

### Rebuild the assay to keep everything consistent
new_assay <- CreateAssayObject(counts = counts_unique)
OMG[[assay_name]] <- new_assay
DefaultAssay(OMG) <- assay_name

# Process
OMG <- NormalizeData(OMG, verbose = FALSE)
OMG <- FindVariableFeatures(OMG, verbose = FALSE)
OMG <- ScaleData(OMG, verbose = FALSE)
OMG <- RunPCA(OMG, npcs = 50, verbose = FALSE)
OMG <- RunUMAP(OMG, dims = 1:50, return.model = T) #, n.components = 3, min.dist = 0.75, 

OMG@meta.data$day <- factor(OMG@meta.data$day, levels = c("E8.0-E8.5","E8.75", "E9.0","E9.25","E9.5","E9.75", "E10.0"))
OMG[["percent.mt"]] <- PercentageFeatureSet(OMG, pattern = "^Mt")

pdf('images/OMG/QC_OMG.pdf')
VlnPlot(OMG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

saveRDS(OMG, 'data/OMG/OMG_E8_E10_merge.rds')

##############################################################################
## UMAP by time point
day_colors <- day_color_plate[c(3:9)]
names(day_colors)[1] <- 'E8.0-E8.5'

pdf('images/OMG/OMG_E8_E10_merge_days.pdf')
DimPlot(OMG, reduction = "umap", group.by='day', raster=T, cols=day_colors)
dev.off()

pdf('images/OMG/OMG_E8_E10_merge_days_split.pdf', width=56)
DimPlot(OMG, reduction = "umap", group.by='day', split.by='day', raster=T, cols=day_colors)
dev.off()

## UMAP by trajectory
pdf('images/OMG/OMG_E8_E10_merge_trajectory_split.pdf', height=9, width=14)
DimPlot(OMG, reduction = "umap", group.by='major_trajectory', label=T, raster=T, cols=major_trajectory_color_plate[levels(factor(OMG@meta.data$major_trajectory))])
dev.off()

## Cell types per day
pdf('images/OMG/OMG_E8_E10_cell_types_per_day.pdf', width=14)
ggplot(OMG@meta.data, aes(x=celltype_update, fill=day)) + geom_bar() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=day_colors)
dev.off()

###############################################################################






