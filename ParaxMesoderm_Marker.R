
# Marker overlap between OMG and TOME cell states (Tbx6+ / Paraxial mesoderm)
#
# Input : tables/scRNAseq/OMG/OMG_filtered_markers.csv  (from 02)
#         tables/scRNAseq/TOME/TOME_filtered_markers.csv (from 03)
#         data/scRNAseq/OMG_filtered.rds and
#         data/scRNAseq/TOME_filtered_20_TOME_E725_E105.rds (shared gene universe)
#
# Steps : 1. Marker overlap between OMG "Mesodermal progenitors (Tbx6+)" and
#            TOME "Paraxial mesoderm A" / "B" (shared-gene counts + percentages,
#            written to CSV)
#         2. Fisher's exact enrichment test for the Tbx6+ / Paraxial A pair,
#            using genes shared across both filtered references as the universe
#         3. All-vs-all OMG x TOME cluster overlap: Fisher's test per pair,
#            BH-adjusted, plotted as a log2 odds-ratio heatmap
#
# Output: tables/scRNAseq/OMG/TOME_parax_meso_{A,B}_..._shared_genes.csv
#         images/scRNAseq/ParaxMeso/OMG_TOME_..._overlap.pdf

##############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

library(ggplot2)

OMG_markers <- read.csv("tables/scRNAseq/OMG/OMG_filtered_markers.csv")
TOME_markers <- read.csv("tables/scRNAseq/TOME/TOME_filtered_markers.csv")

parax_meso_A <- subset(TOME_markers, cluster ==  "Paraxial mesoderm A"  )
tbx6_meso_progenitors <- subset(OMG_markers, cluster == "Mesodermal progenitors (Tbx6+)")

# how many percent are shared between the two lists?
head(parax_meso_A)
shared_genes <- intersect(parax_meso_A$gene, tbx6_meso_progenitors$gene)
percent_shared <- length(shared_genes) / length(tbx6_meso_progenitors$gene) * 100
cat("Counts of shared genes between Paraxial mesoderm A and Tbx6+ Mesodermal progenitors: ", length(shared_genes), "\n")
cat("Percentage of Tbx6+ Mesodermal progenitors that are also in Paraxial mesoderm A: ", round(percent_shared, 2), "% \n")
# export list of shared genes
write.csv(shared_genes, "tables/scRNAseq/OMG/TOME_parax_meso_A_Tbx6_meso_progenitors_shared_genes.csv", row.names = FALSE)

parax_meso_B <- subset(TOME_markers, cluster ==  "Paraxial mesoderm B"  )
tbx6_meso_progenitors <- subset(OMG_markers, cluster == "Mesodermal progenitors (Tbx6+)")

# how many percent are shared between the two lists?
head(parax_meso_B)
shared_genes <- intersect(parax_meso_B$gene, tbx6_meso_progenitors$gene)
percent_shared <- length(shared_genes) / length(tbx6_meso_progenitors$gene) * 100
cat("Counts of shared genes between Paraxial mesoderm B and Tbx6+ Mesodermal progenitors: ", length(shared_genes), "\n")
cat("Percentage of Tbx6+ Mesodermal progenitors that are also in Paraxial mesoderm B: ", round(percent_shared, 2), "% \n")
write.csv(shared_genes, "tables/scRNAseq/OMG/TOME_parax_meso_B_Tbx6_meso_progenitors_shared_genes.csv", row.names = FALSE)

OMG_filtered <- readRDS("data/scRNAseq/OMG_filtered.rds")
integrated_f <- readRDS("data/scRNAseq/TOME_filtered_20_TOME_E725_E105.rds")

universe <- intersect(rownames(OMG_filtered), rownames(integrated_f))

A <- intersect(tbx6_meso_progenitors$gene, universe)   # OMG Tbx6+ markers in universe
B <- intersect(parax_meso_A$gene,          universe)   # TOME Paraxial A markers in universe

shared    <- intersect(A, B)
n_universe <- length(universe)
a <- length(shared)            # in both
b <- length(setdiff(A, B))     # in A only
c <- length(setdiff(B, A))     # in B only
d <- n_universe - a - b - c    # in neither

mat <- matrix(c(a, b, c, d), nrow = 2,
              dimnames = list(InA = c("yes","no"), InB = c("yes","no")))
print(mat)

ft <- fisher.test(mat, alternative = "greater")   # one-sided: enrichment
cat(sprintf("overlap = %d | expected ~%.1f | OR = %.2f | p = %.3g\n",
            a, length(A) * length(B) / n_universe, ft$estimate, ft$p.value))

omg_cl  <- split(OMG_markers$gene,  OMG_markers$cluster)
tome_cl <- split(TOME_markers$gene, TOME_markers$cluster)
universe <- intersect(unique(OMG_markers$gene), unique(TOME_markers$gene))

res <- expand.grid(OMG = names(omg_cl), TOME = names(tome_cl),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
res[c("overlap","OR","p")] <- t(mapply(function(o, t) {
  A <- intersect(omg_cl[[o]],  universe)
  B <- intersect(tome_cl[[t]], universe)
  a <- length(intersect(A, B))
  m <- matrix(c(a, length(A)-a, length(B)-a,
                length(universe)-length(A)-length(B)+a), 2)
  ft <- fisher.test(m, alternative = "greater")
  c(a, unname(ft$estimate), ft$p.value)
}, res$OMG, res$TOME))

res$padj      <- p.adjust(res$p, "BH")          
res$neglog10p <- -log10(res$padj)
res$sig <- ifelse(res$padj < 0.05, "*", "")
res$log2OR <- log2(res$OR)

unique(res$OMG)
unique(res$TOME)

omg_order <- c(
  "Telencephalon",
  "Posterior Forebrain / Diencephalon",
  "Midbrain",
  "Midbrain-hindbrain boundary",
  "Hindbrain",
  "Spinal cord/r7/r8",
  "NMPs and spinal cord progenitors",
  "Mesodermal progenitors (Tbx6+)",
  "Somites / Dermomyotome",
  "Somites / Sclerotome",
  "Notochord",
  "Floorplate and p3 domain",
  "Placodal area",
  "Gut"
)

tome_order <- c(
  "Forebrain",
  "Midbrain",
  "Hindbrain",
  "Spinal cord",
  "Neuromesodermal progenitors",
  "Paraxial mesoderm A",
  "Paraxial mesoderm B",
  "Hematoendothelial progenitors",
  "Endothelium",
  "Gut"
)

res$OMG  <- factor(res$OMG,  levels = omg_order)
res$TOME <- factor(res$TOME, levels = tome_order)

pdf("images/scRNAseq/ParaxMeso/OMG_TOME_parax_meso_A_Tbx6_meso_progenitors_overlap.pdf", width = 5, height = 4)
ggplot(res, aes(x = TOME, y = OMG, fill = log2OR)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(padj < 0.05, "*", "")), size = 4, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log2 odds ratio") +
  coord_fixed() +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = "TOME cluster", y = "OMG cluster")

dev.off()