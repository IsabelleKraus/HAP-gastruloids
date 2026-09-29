# Source Data export: Seurat / Monocle3 objects -> plain tables behind every
# scRNA-seq panel of the manuscript
#
# Input : objects written by 03-06 (data/scRNAseq/)
# Output: source_data/<table>.tsv.gz, read by the scripts in figures/
#
# Run from the repository root.

library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(readr)
library(stringr)
set.seed(42)

source("scripts/OMG_colors.r")                          # celltype_order (OMG numbering)
tome <- new.env(); sys.source("scripts/TOME_colors.r", envir = tome)   # cell_type_ids (TOME numbering)

sd_dir <- "source_data"
dir.create(sd_dir, showWarnings = FALSE)

write_sd <- function(df, name) {
  write.table(df, gzfile(file.path(sd_dir, paste0(name, ".tsv.gz"))),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("%-32s %7d rows", name, nrow(df)))
}

# cell-state numbers as in the figure legends
omg_state_id  <- function(s) as.integer(str_extract(celltype_order[s], "\\d+(?=\\)$)"))
tome_state_id <- function(s) as.integer(str_extract(tome$cell_type_ids[s], "\\d+(?=\\)$)"))

# condition names as in the GEO submission 
recode_cond <- function(x) {
  x <- as.character(x)
  map <- c(Normoxic = "NAP", NAP = "NAP", Hypoxic = "HAP", HAP = "HAP",
           Hypo_XAV = "HAPX", "HAP+X" = "HAPX", HIF1AKO = "HIF1AKO",
           "HAP Hif1a KO" = "HIF1AKO", TLS = "TLS")
  ifelse(x %in% names(map), map[x], x)
}


# 1. TOME (Qiu et al. 2022): Fig 1E, 1F; ED Fig 4A (left)

tome_ref <- readRDS("data/scRNAseq/TOME_filtered_20_TOME_E725_E105.rds")
tome_hap <- readRDS("data/scRNAseq/Asmb_projected_filtered_20_TOME_E725_E105_Fig1F.rds")
tome_tls <- readRDS("data/scRNAseq/TLS_projected_filtered_20_TOME_E725_E105_Fig1F.rds")

emb <- Embeddings(tome_ref, "umap")
cs  <- tome_ref@meta.data[rownames(emb), "cell_type_updated"]
write_sd(data.frame(
  cell_id    = rownames(emb),
  umap_1     = emb[, 1],
  umap_2     = emb[, 2],
  cell_state = cs,
  state_id   = tome_state_id(cs),
  day        = tome_ref@meta.data[rownames(emb), "day"]
), "TOME_reference_UMAP")

tome_query <- function(obj) {
  emb <- Embeddings(obj, "ref.umap")                 # position on the TOME UMAP
  m   <- obj@meta.data[rownames(emb), ]
  data.frame(
    cell_id                    = rownames(emb),
    refumap_1                  = emb[, 1],
    refumap_2                  = emb[, 2],
    condition                  = recode_cond(m$condition),
    cell_state                 = m$predicted.cell_type,
    state_id                   = tome_state_id(m$predicted.cell_type),
    predicted_cell_state_score = m$predicted.cell_type.score,
    predicted_day              = m$predicted.day,
    predicted_day_score        = m$predicted.day.score
  )
}
tq <- rbind(tome_query(tome_tls), tome_query(tome_hap))
print(table(tq$condition))                           # TLS 2438, HAP 1379
write_sd(tq, "TOME_query_UMAP")


# 2. OMG (Qiu et al. 2024): Fig 1G, 1H, 4E, 5E (bars); ED Fig 4A (right), 9E

OMG_filtered <- readRDS("data/scRNAseq/OMG_filtered.rds")
Asmb_OMG     <- readRDS("data/scRNAseq/Asmb_OMG_filtered.rds")   # all conditions, OMG-mapped
TLS_OMG      <- readRDS("data/scRNAseq/TLS_OMG.rds")

emb <- Embeddings(OMG_filtered, "umap")
cs  <- OMG_filtered@meta.data[rownames(emb), "celltype_updated"]
write_sd(data.frame(
  cell_id    = rownames(emb),
  umap_1     = emb[, 1],
  umap_2     = emb[, 2],
  cell_state = cs,
  state_id   = omg_state_id(cs),
  day        = OMG_filtered@meta.data[rownames(emb), "day"]
), "OMG_reference_UMAP")

omg_query <- function(obj) {
  emb <- Embeddings(obj, "ref.umap")                 # position on the OMG UMAP
  m   <- obj@meta.data[rownames(emb), ]
  data.frame(
    cell_id                    = rownames(emb),
    refumap_1                  = emb[, 1],
    refumap_2                  = emb[, 2],
    condition                  = recode_cond(m$condition),
    cell_state                 = m$predicted.celltype_updated,
    state_id                   = omg_state_id(m$predicted.celltype_updated),
    predicted_cell_state_score = m$predicted.celltype_updated.score,
    predicted_day              = m$predicted.day,
    predicted_day_score        = m$predicted.day.score,
    predicted_somite_count     = m$predicted.somite_count
  )
}
oq <- rbind(omg_query(TLS_OMG), omg_query(Asmb_OMG))
print(table(oq$condition, useNA = "ifany"))          # 5 conditions, no NA; HAP 1377
write_sd(oq, "OMG_query_UMAP")
rm(OMG_filtered, TLS_OMG); gc()


# 3. HAP own UMAP + marker expression: Fig 2F-G, 3D, 4E, 5E; ED Fig 5E, 6D

HAP <- readRDS("data/scRNAseq/HAP_OMG.rds")          # HAP UMAP (+ pseudotime from 05)

genes_hap <- c(
  # Fig 2G
  "Emx2", "Pax6", "Otx2", "Otx1", "Dmbx1", "Wnt1", "Egr2", "Gbx2",
  "En1", "En2", "Pax2", "Pax5", "Zic1", "Pax7", "Dbx1", "Shh", "Nkx6-1", "Olig2",
  "Sox10", "Tfap2a",
  # Fig 3D
  "Foxc1", "Foxc2", "Hes7", "Meox1", "Mesp2", "Pax1", "Pax3", "Pax9",
  "T", "Tbx6", "Tbx18", "Uncx",
  # ED Fig 5E
  "Cdx2", "Foxa2", "Gata4", "Gata6", "Hhex", "Sox17"
)
stopifnot(length(setdiff(genes_hap, rownames(HAP))) == 0)

emb  <- Embeddings(HAP, "umap")
cs   <- HAP@meta.data[rownames(emb), "predicted.celltype_updated"]
expr <- as.matrix(GetAssayData(HAP, assay = "RNA", layer = "data")[genes_hap, rownames(emb)])
write_sd(data.frame(
  cell_id    = rownames(emb),
  umap_1     = emb[, 1],
  umap_2     = emb[, 2],
  cell_state = cs,
  state_id   = omg_state_id(cs),
  t(expr),
  check.names = FALSE                               
), "HAP_UMAP_expression")                           


# 4. OMG-predicted somite count per HAP cell: Fig 3C (top)
#    score-weighted sum over somite counts (04: expected_somite)

hap_cells <- colnames(Asmb_OMG)[Asmb_OMG$condition == "HAP"]
somite_scores <- as.matrix(GetAssayData(Asmb_OMG, assay = "prediction.score.somite_count",
                                        layer = "data")[, hap_cells])
somite_scores <- somite_scores[rownames(somite_scores) != "max", , drop = FALSE]
w <- as.numeric(gsub(" somites", "", rownames(somite_scores)))
write_sd(data.frame(
  cell_id                = hap_cells,
  cell_state             = Asmb_OMG@meta.data[hap_cells, "predicted.celltype_updated"],
  predicted_somite_count = colSums(sweep(somite_scores, 1, w, `*`))
), "HAP_predicted_somite_count")


# 5. Score-weighted predicted mouse stage per HAP cell: ED Fig 4B
## Qiu 2024 (OMG) (04, Staging per cell type)
ds <- as.matrix(GetAssayData(Asmb_OMG, assay = "prediction.score.day", layer = "data")[, hap_cells])
ds <- ds[rownames(ds) != "max", , drop = FALSE]
w  <- sapply(rownames(ds), function(x) {             # ranges -> midpoint ("E8.0-E8.5" -> 8.25)
  nums <- as.numeric(regmatches(x, gregexpr("\\d+\\.?\\d*", x))[[1]])
  if (length(nums) == 2) mean(nums) else nums
})
num <- colSums(sweep(ds, 1, w, `*`))
den <- colSums(ds)
omg_stage <- data.frame(
  reference      = "Qiu 2024 (OMG)",
  cell_id        = hap_cells,
  cell_state     = Asmb_OMG@meta.data[hap_cells, "predicted.celltype_updated"],
  expected_stage = num / ifelse(den == 0, 1, den)
)

## Qiu 2022 (TOME) (03, Staging per cell type)
hap_tome <- colnames(tome_hap)[tome_hap$condition == "Hypoxic"]
ds <- as.matrix(GetAssayData(tome_hap, assay = "prediction.score.day", layer = "data")[, hap_tome])
ds <- ds[rownames(ds) != "max", , drop = FALSE]
w  <- parse_number(rownames(ds))                     # "E8.5a" -> 8.5
tome_stage <- data.frame(
  reference      = "Qiu 2022 (TOME)",
  cell_id        = hap_tome,
  cell_state     = tome_hap@meta.data[hap_tome, "predicted.cell_type"],
  expected_stage = colSums(sweep(ds, 1, w, `*`))    
)
write_sd(rbind(tome_stage, omg_stage), "HAP_expected_stage_mouse")
rm(tome_ref, tome_hap, tome_tls, Asmb_OMG); gc()


# 6. NAP and HAP Hif1a KO projected onto HAP UMAP: Fig 4E, 5E-F; ED Fig 6D

Asmb <- readRDS("data/scRNAseq/Asmb_OMG.rds")     

export_projection <- function(cond_internal, code, genes) {
  q <- subset(Asmb, subset = condition == cond_internal)
  q <- NormalizeData(q)
  q <- ScaleData(q, features = rownames(q))
  q <- RunPCA(q, npcs = 50, features = VariableFeatures(q))
  anchors <- FindTransferAnchors(reference = HAP, query = q,
                                 normalization.method = "LogNormalize",
                                 reference.reduction = "pca", dims = 1:30)
  q <- MapQuery(anchorset = anchors, query = q, reference = HAP,
                reference.reduction = "pca", reduction.model = "umap")
  emb  <- Embeddings(q, "ref.umap")                  # position on the HAP UMAP
  cs   <- q@meta.data[rownames(emb), "predicted.celltype_updated"]
  expr <- as.matrix(GetAssayData(q, assay = "RNA", layer = "data")[genes, rownames(emb)])
  data.frame(
    cell_id    = rownames(emb),
    umap_1     = emb[, 1],
    umap_2     = emb[, 2],
    condition  = code,
    cell_state = cs,
    state_id   = omg_state_id(cs),
    t(expr),
    check.names = FALSE
  )
}

write_sd(export_projection("Normoxic", "NAP",
                           c("Emx2", "Pax6", "Otx2", "Dmbx1", "Pax2", "Pax5", "En2", "En1")),
         "NAP_on_HAP_UMAP_expression")
write_sd(export_projection("HIF1AKO", "HIF1AKO",
                           c("Emx2", "Pax6", "Otx1", "Otx2", "Dmbx1", "En1", "En2", "Pax2", "Pax5",
                             "Gata6", "Foxa2", "Sox17")),
         "HIF1AKO_on_HAP_UMAP_expression")


# 7. HAP + XAV own UMAP + markers: Fig 7F-G (same as 06)

HAPX <- subset(Asmb, subset = condition == "Hypo_XAV")
HAPX <- NormalizeData(HAPX)
HAPX <- ScaleData(HAPX, features = rownames(HAPX))
HAPX <- RunPCA(HAPX, npcs = 50, features = VariableFeatures(HAPX))
HAPX <- RunUMAP(HAPX, dims = 1:30, return.model = TRUE)

genes_7G <- c("Emx2", "En1", "Lhx2", "Dmbx1", "Nkx6-1", "Olig2", "Pax3", "Pax7",
              "Sox9", "Pax9", "Foxa2", "Shh", "Tbx18", "Uncx")
stopifnot(length(setdiff(genes_7G, rownames(HAPX))) == 0)

emb  <- Embeddings(HAPX, "umap")
cs   <- HAPX@meta.data[rownames(emb), "predicted.celltype_updated"]
expr <- as.matrix(GetAssayData(HAPX, assay = "RNA", layer = "data")[genes_7G, rownames(emb)])
write_sd(data.frame(
  cell_id    = rownames(emb),
  umap_1     = emb[, 1],
  umap_2     = emb[, 2],
  condition  = "HAPX",
  cell_state = cs,
  state_id   = omg_state_id(cs),
  t(expr),
  check.names = FALSE
), "HAPX_UMAP_expression")

# upper limit of the shared colour scale per gene = max over HAP and HAP + XAV
hap_expr <- as.matrix(GetAssayData(HAP, assay = "RNA", layer = "data")[genes_7G, ])
write_sd(data.frame(gene = genes_7G,
                    vmax = pmax(apply(expr, 1, max), apply(hap_expr, 1, max))),
         "HAP_HAPX_shared_vmax")
rm(Asmb, HAPX); gc()


# 8. Monocle3 trajectory on the HAP UMAP: ED Fig 4E (05, first block)

pt_saved <- HAP$pseudotime                        

Idents(HAP) <- HAP$predicted.celltype_updated
cds_HAP <- as.cell_data_set(HAP)
reducedDims(cds_HAP)$UMAP <- Embeddings(HAP, reduction = "umap")
cds_HAP <- cluster_cells(cds_HAP, reduction_method = "UMAP")
cds_HAP <- learn_graph(cds_HAP, use_partition = TRUE, learn_graph_control = list(ncenter = 200))

g <- principal_graph(cds_HAP)[["UMAP"]]
node_xy <- t(cds_HAP@principal_graph_aux[["UMAP"]]$dp_mst)
colnames(node_xy) <- c("umap_1", "umap_2")

# fork = branch node (degree >= 3) nearest the NMP centroid;
# root = fork neighbour with the highest UMAP 2
branch_nodes <- names(which(igraph::degree(g) >= 3))
nmp_cells    <- colnames(cds_HAP)[colData(cds_HAP)$predicted.celltype_updated ==
                                    "NMPs and spinal cord progenitors"]
nmp_centroid <- colMeans(Embeddings(HAP, "umap")[nmp_cells, ])
d <- sqrt((node_xy[branch_nodes, 1] - nmp_centroid[1])^2 +
          (node_xy[branch_nodes, 2] - nmp_centroid[2])^2)
fork_node <- branch_nodes[which.min(d)]
nbrs      <- names(igraph::neighbors(g, fork_node))
root_node <- nbrs[which.max(node_xy[nbrs, 2])]
message("fork node: ", fork_node, "   root node: ", root_node)

cds_HAP <- order_cells(cds_HAP, root_pr_nodes = root_node)

# same pseudotime as the run saved by 05
if (!is.null(pt_saved)) print(all.equal(unname(pseudotime(cds_HAP)), unname(pt_saved[colnames(cds_HAP)])))

# root direction: fork neighbour through which each cell's path leaves the fork
cv <- as.matrix(cds_HAP@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[
  colnames(cds_HAP), , drop = FALSE])
cell_node      <- setNames(igraph::V(g)$name[as.numeric(cv[, 1])], colnames(cds_HAP))
fork_neighbors <- igraph::V(g)$name[igraph::neighbors(g, fork_node)]

cell_direction <- setNames(rep(NA_character_, length(cell_node)), names(cell_node))
for (cell in names(cell_node)) {
  if (cell_node[cell] == fork_node) { cell_direction[cell] <- "root"; next }
  sp <- igraph::V(g)$name[as.numeric(igraph::shortest_paths(g, from = fork_node, to = cell_node[cell],
                                                            output = "vpath")$vpath[[1]])]
  if (length(sp) >= 2) cell_direction[cell] <- sp[2]
}
# the two largest directions, named by their dominant cell state; the rest -> root
dir_counts <- sort(table(cell_direction[cell_direction %in% fork_neighbors]), decreasing = TRUE)
big2 <- names(dir_counts)[1:2]
dom_ct <- function(nbr) {
  ct <- as.character(colData(cds_HAP)[names(cell_direction)[cell_direction == nbr],
                                      "predicted.celltype_updated"])
  names(sort(table(ct), decreasing = TRUE))[1]
}
somite_nbr <- big2[grepl("Somites", sapply(big2, dom_ct))][1]
brain_nbr  <- setdiff(big2, somite_nbr)[1]
root_direction <- ifelse(cell_direction == brain_nbr,  "direction_Brain",
                  ifelse(cell_direction == somite_nbr, "direction_Somites", "root"))
print(table(root_direction, useNA = "ifany"))

saveRDS(cds_HAP, "data/scRNAseq/HAP_monocle_cds_EDFig4E.rds")

cs <- as.character(colData(cds_HAP)$predicted.celltype_updated)
pt <- pseudotime(cds_HAP)
write_sd(data.frame(
  cell_id        = colnames(cds_HAP),
  umap_1         = reducedDims(cds_HAP)$UMAP[, 1],
  umap_2         = reducedDims(cds_HAP)$UMAP[, 2],
  cell_state     = cs,
  state_id       = omg_state_id(cs),
  pseudotime     = ifelse(is.finite(pt), pt, NA),
  root_direction = unname(root_direction)
), "HAP_trajectory_cells")

el <- igraph::as_edgelist(g)                         
write_sd(data.frame(from = el[, 1], to = el[, 2],
                    x    = node_xy[el[, 1], 1], y    = node_xy[el[, 1], 2],
                    xend = node_xy[el[, 2], 1], yend = node_xy[el[, 2], 2]),
         "HAP_trajectory_graph")


writeLines(capture.output(sessionInfo()), file.path(sd_dir, "sessionInfo_export.txt"))