# HAP-gastruloids — single-cell RNA-seq analysis

Analysis code for the single-cell RNA-seq study of antero-posterior (AP)
gastruloids grown under hypoxia (HAP-gastruloids). ESC-derived structures were
profiled across four conditions — normoxic (NAP), hypoxic (HAP), hypoxic + XAV939
(HAPX), and HIF1A-KO and annotated against mouse (Qiu et al. 2022, Qiu et al. 2024) and human
(Xu et al. 2023, Zeng et al. 2023) developmental references.

> The developing embryo is guided by continuously changing signals from its
> microenvironment, among these, restricted oxygen (hypoxia) is a critical
> regulator of cell-type diversification. This work builds a modular ESC-based
> head-to-tail model of mouse embryogenesis under hypoxia. The resulting
> HAP-gastruloids form stage-appropriate anterior neural tissue (fore-/midbrain,
> including the midbrain–hindbrain boundary) in synchrony with posterior tissues
> (spinal cord, somites, gut endoderm derivatives), and genetic, environmental,
> and pharmacological perturbations show that timed hypoxia is essential for
> forebrain identity and proper neural patterning.

---

## Data availability

Processed count matrices are deposited at **GEO: GSE337584** (four Parse
Biosciences conditions + the TLS 10x dataset).

The GEO matrices are condition-prefixed files, three per condition:

```
NAP_count_matrix.mtx.gz     NAP_all_genes.csv.gz     NAP_cell_metadata.csv.gz
HAP_count_matrix.mtx.gz     HAP_all_genes.csv.gz     HAP_cell_metadata.csv.gz
HAPX_count_matrix.mtx.gz    HAPX_all_genes.csv.gz    HAPX_cell_metadata.csv.gz
HIF1AKO_count_matrix.mtx.gz HIF1AKO_all_genes.csv.gz HIF1AKO_cell_metadata.csv.gz
```

**Condition-name mapping (GEO → internal script names):**

| GEO | Script | Meaning |
|-----|--------|---------|
| NAP | Normoxic | normoxic |
| HAP | Hypoxic | hypoxic |
| HAPX | Hypo_XAV | hypoxic + XAV939 |
| HIF1AKO | HIF1AKO | HIF1A knockout |

To run `04` on the GEO files, place each condition's three files under
`data/scRNAseq/filtered_matrices/output_combined/<condition>/DGE_filtered/`,
renamed to `count_matrix.mtx.gz` / `cell_metadata.csv.gz` / `all_genes.csv.gz`
(or edit the `base` path and `read_parse_sample()` at the top of `04` to read the
flat GEO names directly). Public reference atlases (OMG/TOME/Xu/Zeng) are
downloaded from their original sources. See the header comments in `01`, `03`,
and `07` for the download URLs and accessions.

---

## Repository layout

All paths in the scripts are **relative to the repository root**. Set the
working directory to the repo root before running.
Large inputs and generated outputs are not tracked (see
`.gitignore`); create the `data/`, `images/`, and `tables/` folders as needed.

```
01_build_OMG_reference.R      Build the OMG reference (Qiu et al. 2024)
02_annotate_OMG_reference.R   Markers, subclustering
03_build_TOME_reference.R     Build the TOME reference (Qiu et al. 2022), subcluster
                              neural states, stage the gastruloids
04_process_gastruloids_OMG.R  Main: load Parse matrices, QC + scDblFinder,
                              OMG label transfer, composition / staging / markers
05_pseudotime.R               Monocle3 pseudotime trajectory (HAP condition)
06_marker_expression.R        Per-condition marker FeaturePlots and
                              cross-condition UMAP projections
07_human_references.R         Cross-species staging vs. human atlases
                              (Xu 2023, Zeng 2023)
ParaxMesoderm_Marker.R        Supplementary: OMG and TOME marker-overlap analysis
scripts/OMG_colors.r          Colour palettes (OMG cell states / stages)
scripts/TOME_colors.r         Colour palettes (TOME cell states / stages)
app/                          HAPView — interactive Shiny explorer (see below)
```

**Note on ordering:** `04` writes an intermediate QC'd object
(`data/scRNAseq/Asmb.rds`) early, which `03` reads for its TOME staging step;
In practice the dependency is: `04` (QC) → `03`

---

## Requirements

R 4.4.1. Core packages:

```r
install.packages(c("Seurat", "Matrix", "dplyr", "tibble", "tidyr", "ggplot2",
                   "stringr", "patchwork", "circlize", "ComplexHeatmap",
                   "ggrastr", "ggrepel", "RColorBrewer"))
# trajectory / integration / references:
#   monocle3, SeuratWrappers, scDblFinder, DoubletFinder, homologene,
#   org.Mm.eg.db, AnnotationDbi
```

Key versions used: Seurat 5.1.0, Monocle3, scDblFinder 1.18.0. Parse data were
processed upstream with Trailmaker (split-pipe v1.7.3); TLS data with Cell Ranger.

---

## HAPView — interactive Shiny app

`app/` contains **HAPView**, a Shiny app for exploring any Seurat `.rds` object.

| Tab | Contents |
|-----|----------|
| UMAP / PCA | Reductions coloured by cell type and by gene expression; composition bar chart |
| Heatmap | Average-expression heatmap (auto `FindAllMarkers()` or custom gene list) |
| Gene Expression | Per-cell-type violin plots + multi-gene dot plots |
| Annotations | Cell-type counts, pie chart, metadata browser |
| Dataset Info | Seurat object summary |

**Input:** Seurat v4/v5 `.rds` with at least one reduction (`umap`/`pca`) and a
cell-type column (auto-detects names containing *celltype* / *cluster* / *ident*
/ *annotation*; for this dataset, `predicted.celltype_updated`).

**Dependencies** (installed automatically on first launch if missing):
`shiny`, `shinydashboard`, `shinyWidgets`, `DT`, `Seurat`, `ggplot2`, `plotly`,
`viridis`, `pheatmap`, `RColorBrewer`, `dplyr`, `tibble`, `scales`, `ggrepel`,
`patchwork`.

**Run it** (the app file is `app/app.R`):

```r
# RStudio: open app/app.R and click "Run App", or:
shiny::runApp("app")
```

```bash
# On a server:
Rscript -e "shiny::runApp('app', port = 3838, host = '0.0.0.0')"
# then open http://<your-server>:3838

# Optionally auto-load an object on startup:
TASOSOID_PRELOAD="data/HAP_OMG.rds" \
  Rscript -e "shiny::runApp('app', port = 3838, host = '0.0.0.0')"
```

Performance notes: the metadata table previews the first 500 cells (full data used for counts/pie); `FindAllMarkers()` on a large object can take 1–5 min.

---

## Citation

If you use this code, please cite the associated publication *(add reference /
DOI on release)* and the GEO accession above.
