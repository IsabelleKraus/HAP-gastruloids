
# Pseudotime / trajectory analysis of the HAP condition (Monocle3)
#
# Input : data/scRNAseq/Asmb_OMG.rds and data/scRNAseq/HAP_OMG.rds
#         (HAP subset with OMG cell-state annotation, from 04)
#
# Steps : 1. Convert HAP to a Monocle3 cell_data_set, carrying over the Seurat UMAP 
#         2. cluster_cells() + learn_graph(use_partition = TRUE, ncenter = 200)
#         3. Programmatic root: fork = principal-graph node of degree >=3
#            nearest the NMP / spinal-cord-progenitor centroid; root = fork
#            neighbour with highest UMAP-2 -> order_cells()
#         4. Save HAP with a pseudotime column
#         5. Pseudotime heatmaps (binned + per-cell)
#         6. Split cells into neural / mesodermal branches by shortest path to
#            the fork; per-branch heatmaps + marker line plots
#
# Output: data/scRNAseq/HAP_OMG.rds (with pseudotime),
#         figures in images/scRNAseq/Monocle/

###############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

library(monocle3)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(magrittr)
library(ggplot2)
library(Seurat)
library(future)
library(dplyr)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(tidyr)
set.seed(42)
source("scripts/OMG_colors.r")
Asmb <- readRDS("data/scRNAseq/Asmb_OMG.rds")
HAP <- readRDS("data/scRNAseq/HAP_OMG.rds")

Idents(HAP) <- HAP$predicted.celltype_updated
cds_HAP <- as.cell_data_set(HAP)
umap_coords <- Embeddings(HAP, reduction = "umap")
# Add UMAP coordinates to Monocle object
reducedDims(cds_HAP)$UMAP <- umap_coords

cds_HAP <- cluster_cells(cds_HAP, reduction_method = "UMAP")
cds_HAP <- learn_graph(cds_HAP, use_partition = TRUE, learn_graph_control = list(
                           ncenter = 200))

# Plot the trajectory graph
  pdf("images/scRNAseq/Monocle/HAP_trajectory_states.pdf", width=10, height=7)
  p1 <- plot_cells(cds_HAP, cell_size=3,label_groups_by_cluster=FALSE, label_cell_groups = FALSE,trajectory_graph_color="black",rasterize = FALSE,color_cells_by="predicted.celltype_updated",
            label_leaves=FALSE, 
            label_branch_points=FALSE) + 
    scale_color_manual(values=celltype_updated_colors) + guides(color = guide_legend(override.aes = list(stroke = 0, color=NA))) + coord_fixed(ratio=1) 
  p1$layers[[1]]$aes_params$colour <- 'transparent'
  p1 <- p1 + theme_classic() +
    theme(
      panel.border     = element_blank(), 
      axis.line        = element_blank(), 
      axis.text.x      = element_blank(), 
      axis.text.y      = element_blank(), 
      axis.ticks       = element_blank(), 
      axis.title.x     = element_blank(), 
      axis.title.y     = element_blank(), 
      panel.background = element_blank()  
    )
    print(p1)
  dev.off()


g <- principal_graph(cds_HAP)[["UMAP"]]
node_xy <- t(cds_HAP@principal_graph_aux[["UMAP"]]$dp_mst)
colnames(node_xy) <- c("umap_1", "umap_2")
if (nrow(node_xy) == 2) node_xy <- t(node_xy)   # guard orientation

# fork node
branch_nodes <- names(which(igraph::degree(g) >= 3))
nmp_cells <- colnames(cds_HAP)[colData(cds_HAP)$predicted.celltype_updated ==
                               "NMPs and spinal cord progenitors"]
nmp_centroid <- colMeans(Embeddings(HAP, "umap")[nmp_cells, ])
d <- sqrt((node_xy[branch_nodes,1] - nmp_centroid[1])^2 +
          (node_xy[branch_nodes,2] - nmp_centroid[2])^2)
fork_node <- branch_nodes[which.min(d)]

# neighbors of the fork, pick the one with the highest umap_2
nbrs <- names(igraph::neighbors(g, fork_node))
root_node <- nbrs[which.max(node_xy[nbrs, 2])]
root_node

cds_HAP <- order_cells(cds_HAP, root_pr_nodes = root_node)
umap_df <- as.data.frame(reducedDims(cds_HAP)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$celltype <- colData(cds_HAP)$predicted.celltype_updated
centroids <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype, data = umap_df, FUN = median)

pt <- pseudotime(cds_HAP)
stopifnot(all(names(pt) == colnames(HAP)))  
HAP$pseudotime <- pt
saveRDS(HAP, "data/scRNAseq/HAP_OMG.rds")

# Plot pseudotime and save
pdf("images/scRNAseq/Monocle/HAP_trajectory_pseudotime.pdf", width=10, height=7)
p2 <- plot_cells(cds_HAP,
                 color_cells_by         = "pseudotime",
                 label_leaves           = FALSE,
                 label_roots            = TRUE,
                 label_groups_by_cluster = FALSE,
                 cell_size = 3,
                 label_cell_groups      = FALSE,
                 label_branch_points    = FALSE,
                 rasterize = FALSE,
                 trajectory_graph_color = "black") +
  geom_text_repel(data = centroids,
                  aes(x = UMAP_1, y = UMAP_2, label = celltype),
                  size = 6, max.overlaps = Inf,
                  bg.color = "white", bg.r = 0.15,
                  inherit.aes = FALSE) + 
  coord_fixed(ratio=1)+ 
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
  p2$layers <- p2$layers[-1]
p2 <- p2 + theme_classic()+
  theme(
    panel.border     = element_blank(), 
    axis.line        = element_blank(), 
    axis.text.x      = element_blank(), 
    axis.text.y      = element_blank(), 
    axis.ticks       = element_blank(), 
    axis.title.x     = element_blank(), 
    axis.title.y     = element_blank(), 
    panel.background = element_blank()  
  )
print(p2)
dev.off()

#########################################################################

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
markers_pt <- markers[markers %in% rownames(HAP)]

# remove cells without pseudotime 
pt <- pseudotime(cds_HAP)
cells_finite  <- names(pt[is.finite(pt)])
pt_ordered    <- sort(pt[cells_finite])
cells_ordered <- names(pt_ordered)

# bin cells into pseudotime windows and average expression per bin
n_bins   <- 100
bin_idx  <- cut(seq_along(cells_ordered), breaks = n_bins, labels = FALSE)
expr_mat <- as.matrix(GetAssayData(HAP, layer = "data")[markers_pt, cells_ordered])

expr_binned <- do.call(cbind, lapply(seq_len(n_bins), function(b) {
  idx <- which(bin_idx == b)
  if (length(idx) == 1) expr_mat[, idx] else rowMeans(expr_mat[, idx])
}))
colnames(expr_binned) <- seq_len(n_bins)
expr_scaled <- t(scale(t(expr_binned)))

# majority cell type per bin 
celltypes_ordered <- colData(cds_HAP)[cells_ordered, "predicted.celltype_updated"]
bin_celltype <- sapply(seq_len(n_bins), function(b) {
  idx <- which(bin_idx == b)
  ct  <- as.character(celltypes_ordered[idx])
  names(sort(table(ct), decreasing = TRUE))[1]
})

keep <- rowSums(is.finite(expr_scaled)) == ncol(expr_scaled)
message("dropping flat markers: ", paste(rownames(expr_scaled)[!keep], collapse = ", "))
expr_scaled <- expr_scaled[keep, ]

ha_top <- HeatmapAnnotation(
  celltype = bin_celltype,
  col      = list(celltype = celltype_updated_colors),
  annotation_name_side = "left",
  show_legend = TRUE
)

ht <- Heatmap(
  expr_scaled,
  name              = "z-score",
  top_annotation    = ha_top,
  cluster_columns   = FALSE,
  cluster_rows      = FALSE,
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 8),
  col               = colorRamp2(c(min(expr_scaled), 0, 2, max(expr_scaled)), c("#2166AC", "white", "#B2182B", "#B2182B")),
  use_raster        = TRUE,
  column_title      = "Pseudotime →"
)

pdf("images/scRNAseq/Monocle/HAP_pseudotime_heatmap.pdf", width = 10, height = 10)
draw(ht)
dev.off()

# Pseudotime heatmap (per-cell)
expr_scaled_cells <- t(scale(t(expr_mat)))

ha_top_cells <- HeatmapAnnotation(
  celltype = as.character(celltypes_ordered),
  col      = list(celltype = celltype_updated_colors),
  annotation_name_side = "left",
  show_legend = TRUE
)

ht_cells <- Heatmap(
  expr_scaled_cells,
  name              = "z-score",
  top_annotation    = ha_top_cells,
  cluster_columns   = FALSE,
  cluster_rows      = FALSE,
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 8),
  col               = colorRamp2(c(min(expr_scaled_cells, na.rm=TRUE), 0, 2, max(expr_scaled_cells, na.rm=TRUE)),
                                 c("#2166AC", "white", "#B2182B", "#B2182B")),
  use_raster        = TRUE,
  column_title      = "Pseudotime →"
)

pdf("images/scRNAseq/Monocle/HAP_pseudotime_heatmap_cells.pdf", width = 10, height = 8)
draw(ht_cells)
dev.off()


# Split cells by branch from NMPs 
selected_root <- fork_node 
# Split cells by directions from the already selected NMP root 
pg <- principal_graph(cds_HAP)[["UMAP"]]
igraph::degree(pg, v = selected_root)

# closest principal node per cell
cv <- as.matrix(
  cds_HAP@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[
    colnames(cds_HAP), , drop = FALSE]
)
cell_node <- igraph::V(pg)$name[as.numeric(cv[, 1])]
names(cell_node) <- colnames(cds_HAP)

# the fork-neighbor each cell's path leaves through
fork_neighbors <- igraph::V(pg)$name[igraph::neighbors(pg, selected_root)]
cell_direction <- rep(NA_character_, length(cell_node))
names(cell_direction) <- names(cell_node)

for (cell in names(cell_node)) {
  this_node <- cell_node[cell]
  if (this_node == selected_root) { cell_direction[cell] <- "root"; next }
  sp <- igraph::shortest_paths(pg, from = selected_root, to = this_node,
                               output = "vpath")$vpath[[1]]
  sp_names <- igraph::V(pg)$name[as.numeric(sp)]
  if (length(sp_names) >= 2) cell_direction[cell] <- sp_names[2]  # first step
}

# rank paths by size, keep the 2 biggest
dir_counts <- sort(table(cell_direction[cell_direction %in% fork_neighbors]),
                   decreasing = TRUE)
print(dir_counts)
big2 <- names(dir_counts)[1:2]

# dominant cell type of 2 branches
dom_ct <- function(nbr) {
  cells <- names(cell_direction)[cell_direction == nbr]
  ct <- as.character(colData(cds_HAP)[cells, "predicted.celltype_updated"])
  names(sort(table(ct), decreasing = TRUE))[1]
}
print(sapply(big2, dom_ct))  

somite_nbr <- big2[grepl("Somites",
                         sapply(big2, dom_ct))][1]
brain_nbr  <- setdiff(big2, somite_nbr)[1]
final <- ifelse(cell_direction == brain_nbr,  "direction_neural",
         ifelse(cell_direction == somite_nbr, "direction_mesodermal", "root"))

colData(cds_HAP)$root_direction <- factor(
  final, levels = c("direction_neural", "direction_mesodermal", "root"))
colData(cds_HAP)$branch <- as.character(colData(cds_HAP)$root_direction)
table(colData(cds_HAP)$root_direction, useNA = "ifany")

cols <- c("root" = "#9B1C31", "direction_neural" = "#EAA448", "direction_mesodermal" = "#92B9BD")

# Extract pseudotime values
cell_p_time <- pseudotime(cds_HAP)
# Sort cells from highest to lowest pseudotime
cds_sorted_early <- cds_HAP[, order(cell_p_time, decreasing = TRUE)]

pdf("images/scRNAseq/Monocle/HAP_trajectory_root_directions.pdf", width=8, height=6)
p <- plot_cells(
  cds_sorted_early,
  color_cells_by = "root_direction",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  
  label_branch_points = FALSE,
  label_roots = FALSE,
  cell_size = 3
)+
scale_color_manual(values = cols) +
  theme_classic()+
  theme(
    panel.border     = element_blank(), 
    axis.line        = element_blank(), 
    axis.text.x      = element_blank(), 
    axis.text.y      = element_blank(), 
    axis.ticks       = element_blank(), 
    axis.title.x     = element_blank(), 
    axis.title.y     = element_blank(), 
    panel.background = element_blank()  
  )
p$layers[[1]]$aes_params$colour <- 'transparent'
print(p)
dev.off()

# Use root-based directions as branches
branch_label <- as.character(colData(cds_HAP)$root_direction)
colData(cds_HAP)$branch <- branch_label
table(colData(cds_HAP)$branch, useNA = "ifany")

dominant_ct <- function(cells) {
  ct <- as.character(colData(cds_HAP)[cells, "predicted.celltype_updated"])
  names(sort(table(ct), decreasing = TRUE))[1]
}

cat(
  "Branch Neural dominant type:",
  dominant_ct(rownames(colData(cds_HAP))[colData(cds_HAP)$branch == "direction_neural"]),
  "\n"
)
# Branch Neural dominant type: Midbrain-hindbrain boundary 
cat(
  "Branch Mesodermal dominant type:",
  dominant_ct(rownames(colData(cds_HAP))[colData(cds_HAP)$branch == "direction_mesodermal"]),
  "\n"
)
# Branch Mesodermal dominant type: Somites / Dermomyotome 

make_branch_plots <- function(
    branch_name,
    seurat_obj,
    cds_obj,
    markers_use,
    celltype_colors,
    out_dir,
    include_root = TRUE
) {
  
  # Select cells from this branch
  cells_branch <- rownames(colData(cds_obj))[
    colData(cds_obj)$branch == branch_name
  ]
  
  # Optionally include root cells at the beginning of both branches
  if (include_root) {
    cells_root <- rownames(colData(cds_obj))[
      colData(cds_obj)$branch == "root"
    ]
    cells <- unique(c(cells_root, cells_branch))
  } else {
    cells <- cells_branch
  }
  
  pt_b <- pseudotime(cds_obj)[cells]
  pt_b <- pt_b[is.finite(pt_b)]
  cells_ord <- names(sort(pt_b))
  
  if (length(cells_ord) < 10) {
    warning("Too few cells for branch: ", branch_name)
    return(NULL)
  }
  
  markers_b <- markers_use[markers_use %in% rownames(seurat_obj)]
  
  # Number of bins should not exceed number of cells
  n_b <- min(100, length(cells_ord))
  
  bin_b <- cut(
    seq_along(cells_ord),
    breaks = n_b,
    labels = FALSE
  )
  
  emat <- as.matrix(
    GetAssayData(seurat_obj, layer = "data")[markers_b, cells_ord, drop = FALSE]
  )
  
  ebinned <- do.call(cbind, lapply(seq_len(n_b), function(b) {
    idx <- which(bin_b == b)
    if (length(idx) == 1) {
      emat[, idx]
    } else {
      rowMeans(emat[, idx, drop = FALSE])
    }
  }))
  
  colnames(ebinned) <- seq_len(n_b)
  
  # Scale expression per gene
  escalated <- t(scale(t(ebinned)))
  # Drop genes with zero variance 
  keep <- rowSums(is.na(escalated)) == 0

  escalated <- escalated[keep, , drop = FALSE]
  ebinned   <- ebinned[keep, , drop = FALSE]
  
  if (nrow(ebinned) == 0) {
    warning("No variable genes left after scaling for branch: ", branch_name)
    return(NULL)
  }
  
  # Majority cell type per bin
  cts_ord <- colData(cds_obj)[cells_ord, "predicted.celltype_updated"]
  
  bin_ct_b <- sapply(seq_len(n_b), function(b) {
    ct <- as.character(cts_ord[which(bin_b == b)])
    names(sort(table(ct), decreasing = TRUE))[1]
  })
  
  # Keep only colors present in this branch
  celltype_colors_use <- celltype_colors[
    names(celltype_colors) %in% unique(bin_ct_b)
  ]

  ha <- HeatmapAnnotation(
    celltype = bin_ct_b,
    col = list(celltype = celltype_colors),
    annotation_name_side = "left",
    show_legend = TRUE
  )
  
  ht_b <- Heatmap(
    escalated,
    name = "z-score",
    top_annotation = ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    col = colorRamp2(
      c(min(escalated, na.rm = TRUE), 0, 5, max(escalated, na.rm = TRUE)),
      c("#2166AC", "white", "#B2182B", "#B2182B")
    ),
    use_raster = TRUE,
    column_title = paste0(branch_name, " — Pseudotime →")
  )
  
  pdf(
    file.path(out_dir, paste0("HAP_", branch_name, "_heatmap.pdf")),
    width = 10,
    height = 10
  )
  draw(ht_b)
  dev.off()
  
  # Line plots: normalize per gene to 0–1
  enorm <- t(apply(ebinned, 1, function(x) {
    r <- range(x, na.rm = TRUE)
    if (diff(r) == 0) {
      rep(0, length(x))
    } else {
      (x - r[1]) / diff(r)
    }
  }))
  
  colnames(enorm) <- colnames(ebinned)
  
gene_label <- rownames(enorm)
gene_id_order <- make.unique(gene_label, sep = "_dup")
rownames(enorm) <- gene_id_order

elong <- as.data.frame(enorm) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(
    -gene_id,
    names_to = "bin",
    values_to = "expression"
  ) %>%
  mutate(
    bin = as.integer(bin),
    gene_id = factor(gene_id, levels = gene_id_order)
  )

gene_labels_named <- setNames(gene_label, gene_id_order)

  ct_df_b <- data.frame(
    bin = seq_len(n_b),
    celltype = factor(bin_ct_b, levels = names(celltype_colors))
  )
  ct_long_b <- ct_df_b[
  rep(seq_len(n_b), length(gene_id_order)),
]

ct_long_b$gene_id <- factor(
  rep(gene_id_order, each = n_b),
  levels = gene_id_order
)
  
  p_lp <- ggplot(elong, aes(x = bin, y = expression)) +
    geom_tile(
      data = ct_long_b,
      aes(x = bin, y = -0.12, fill = celltype),
      height = 0.1,
      inherit.aes = FALSE
    ) +
    geom_line(linewidth = 0.5) +
    scale_fill_manual(
      values = celltype_colors,
      name = "Cell type",
      na.value = "grey80"
    ) +
    facet_wrap(
      ~gene_id,
      ncol = 6,
      labeller = labeller(gene_id = gene_labels_named)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(-0.2, 1.05),
      breaks = c(0, 0.5, 1)
    ) +
    labs(
      x = "Pseudotime ->",
      y = "Norm. expression",
      title = paste0(branch_name, " branch")
    ) +
    theme_classic(base_size = 8) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 7),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.3, "lines")
    )
  
  pdf(
    file.path(out_dir, paste0("HAP_", branch_name, "_lineplots.pdf")),
    width = 12,
    height = 12
  )
  print(p_lp)
  dev.off()
  
  return(list(
    heatmap = ht_b,
    lineplot = p_lp,
    cells = cells_ord,
    bin_celltype = bin_ct_b
  ))
}


out_dir <- "images/scRNAseq/Monocle"

res_A <- make_branch_plots(
  branch_name = "direction_neural",
  seurat_obj = HAP,
  cds_obj = cds_HAP,
  markers_use = markers_pt,
  celltype_colors = celltype_updated_colors,
  out_dir = out_dir,
  include_root = TRUE
)

res_B <- make_branch_plots(
  branch_name = "direction_mesodermal",
  seurat_obj = HAP,
  cds_obj = cds_HAP,
  markers_use = markers_pt,
  celltype_colors = celltype_updated_colors,
  out_dir = out_dir,
  include_root = TRUE
)




















































