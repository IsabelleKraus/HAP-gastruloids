
#####################################################################################
# HAPView — Shiny app for interactive exploration of Seurat scRNA-seq objects
# Tabs: UMAP/PCA, Heatmap (FindAllMarkers or custom genes), Gene Expression
#       (violin / dot plots), Annotations, Dataset Info
#
# Input : any Seurat v4/v5 .rds with a reduction (umap/pca) and a cell-type
#         column (auto-detects names containing celltype/cluster/ident/annotation)
#
# Run   : open in RStudio and click "Run App", or
#           shiny::runApp("app.R", port = 3838, host = "0.0.0.0")
#         Optionally auto-load an object on startup:
#           TASOSOID_PRELOAD=/path/to/object.rds Rscript -e 'shiny::runApp(...)'
#####################################################################################

#  HAPView — Single-Cell RNA-seq Visualization Shiny App
#  Supports: Seurat v4/v5 RDS objects
#  Plots: UMAP, PCA, Heatmap, Cell-type annotation summary

# Dependencies 
required_pkgs <- c(
  "shiny", "shinydashboard", "shinyWidgets", "DT",
  "Seurat", "ggplot2", "plotly", "viridis",
  "pheatmap", "RColorBrewer", "dplyr", "tibble",
  "scales", "ggrepel", "patchwork"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(DT)
  library(Seurat)
  library(ggplot2)
  library(plotly)
  library(viridis)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(scales)
  library(ggrepel)
})


#  Cell type color palettes 
# Primary palette for predicted.celltype_updated (and general use)
celltype_updated_colors <- c(
  "Anterior floor plate" = "#a14a57", "Floorplate and p3 domain" = "#a0522b",
  "Anterior roof plate" = "#FFA69E", "Posterior roof plate" = "#ffb6a0",
  "Spinal cord ventral progenitors" = "#e97451", "Spinal cord motor neurons" = "#ffaf87",
  "Spinal cord dorsal progenitors" = "#e3b778", "Spinal cord/r7/r8" = "#ffdcad",
  "NMPs and spinal cord progenitors" = "#9b1c31", "Dorsal telencephalon" = "#ff7f50",
  "Telencephalon" = "#c50061", "Posterior Forebrain / Diencephalon" = "red",
  "Midbrain" = "#c1996b", "Hindbrain" = "#d94e1f",
  "Midbrain-hindbrain boundary" = "#ffae42", "Hypothalamus" = "#CC3131",
  "Hypothalamus (Sim1+)" = "#D34A4A", "Anterior Forebrain" = "#ff4500",
  "Cranial motor neurons" = "#dc143c", "Neural crest (PNS neurons)" = "#732507",
  "Neural crest (PNS glia)" = "#8b4513",
  "Olfactory epithelial cells" = "#ff1493", "Otic epithelial cells" = "#ff69b4",
  "Pre-epidermal keratinocytes" = "#ffb6c1",
  "Somites / Dermomyotome" = "#97e2ff", "Somites / Sclerotome" = "#156888",
  "Mesodermal progenitors (Tbx6+)" = "#6689a1",
  "Anterior intermediate mesoderm" = "#4f4a83",
  "Lateral plate and intermediate mesoderm" = "#1d3557",
  "First heart field" = "#c600cb", "Second heart field" = "#800080",
  "Endocardial cells" = "#E3B5E3", "Arterial endothelial cells" = "#da70d6",
  "Endothelium" = "#d1f1d2", "Pericytes" = "#90ee90",
  "Hematoendothelial progenitors" = "#6c8e23", "Primitive erythroid cells" = "#556b2f",
  "Facial mesenchyme" = "#4b0082", "FM_1" = "#6e5586", "FM_0" = "#b084c0",
  "FM_2" = "#4b0082", "Chondrocytes (Atp1a2+)" = "#8a2be2",
  "Limb mesenchyme progenitors" = "#9370db",
  "Gut" = "#f6f201", "Ciliated nodal cells" = "#ffd700",
  "Primordial germ cells" = "#ff8dfc",
  "Notochord" = "#1242fd", "Olfactory pit cells" = "#1e90ff",
  "Granular keratinocytes" = "#87cefa", "Pancreatic acinar cells" = "#00bfff",
  "Pituitary/Pineal gland progenitors" = "#4682b4", "Apical ectodermal ridge" = "#5f9ea0",
  "Amniotic ectoderm" = "#ffa07a", "Placodal area" = "#DDA0DD",
  "Posterior intermediate mesoderm" = "#7B68EE",
  "Otic sensory neurons" = "#FF69B4", "Sympathetic neurons" = "#FF47A9",
  "Enteric neurons" = "#A30057", "Olfactory sensory neurons" = "#FF26FF",
  "Muscle progenitor cells" = "#cd5c5c", "Muscle progenitor cells (Prdm1+)" = "#b22222",
  "Myotubes" = "#8b0000", "Hepatocytes" = "#daa520",
  "Lung progenitor cells" = "#ffdead", "Pancreatic islets" = "#20B2AA",
  "Suprachiasmatic nucleus" = "#008080", "Megakaryocytes" = "#2E8B57",
  "Hematopoietic stem cells (Cd34+)" = "#006400",
  "Definitive early erythroblasts (CD36-)" = "#228B22",
  "Border-associated macrophages" = "#32CD32",
  "Border-associated macrophages (Ms4a8a+)" = "#3CB371",
  "GABAergic neurons" = "#783B3B", "GABAergic cortical interneurons" = "#5C2A2A",
  "Glutamatergic neurons" = "#DB6B6B", "Neural progenitor cells (Neurod1+)" = "#E38C8C",
  "Choroid plexus" = "#6C350F", "Brain capillary endothelial cells" = "#AF5931",
  "Liver sinusoidal endothelial cells" = "#FFA500",
  "Extraembryonic visceral endoderm" = "#f6c100",
  "Midgut/Hindgut epithelial cells" = "#ffecb3"
)

# Palette for predicted.day
day_colors <- c(
  "E8.0-E8.5" = "#78c1e2",
  "E8.75"     = "#2c7bb6",
  "E9.0"      = "#75b38d",
  "E9.25"     = "#2e5d37",
  "E9.5"      = "#c786cd",
  "E9.75"     = "#8f3b96"
)

# Reproducible random palette generator (same labels always get same colors)
random_palette <- function(labels, seed_offset = 0) {
  lvls <- unique(as.character(labels))
  # Deterministic seed derived from label set so colors are stable across sessions
  set.seed(sum(utf8ToInt(paste(sort(lvls), collapse = ""))) + seed_offset)
  hues <- sample(360, length(lvls), replace = TRUE)
  cols <- hcl(h = hues, c = 65, l = 65)
  setNames(cols, lvls)
}

# Helper: resolve colors for a character vector of labels.
# Automatically picks the best palette:
#   1. If labels match day palette (e.g. "E8.75") -> day_colors
#   2. If labels match celltype palette -> celltype_updated_colors
#   3. Otherwise -> reproducible random colors
# Any unmatched individual labels within a chosen palette get random-filled.
resolve_colors <- function(labels, palette = NULL) {
  lvls <- unique(as.character(labels))

  # Auto-pick palette when not explicitly provided
  if (is.null(palette)) {
    day_hits <- sum(lvls %in% names(day_colors))
    ct_hits  <- sum(lvls %in% names(celltype_updated_colors))
    palette  <- if (day_hits >= ct_hits && day_hits > 0) day_colors
                else if (ct_hits > 0)                    celltype_updated_colors
                else                                     NULL   # no palette match -> all random
  }

  if (is.null(palette)) {
    return(random_palette(lvls))
  }

  # Match what we can from the palette; random-fill anything unknown
  cols <- setNames(rep(NA_character_, length(lvls)), lvls)
  known <- lvls %in% names(palette)
  cols[known] <- palette[lvls[known]]
  if (any(!known)) {
    fill_cols <- random_palette(lvls[!known])
    cols[!known] <- fill_cols[lvls[!known]]
  }
  cols
}

# Pre-load config 
# Set env var TASOSOID_PRELOAD to an RDS path to auto-load on startup
DEFAULT_PRELOAD <- Sys.getenv("TASOSOID_PRELOAD", unset = "")
PRELOAD_MODE    <- nzchar(DEFAULT_PRELOAD) && file.exists(DEFAULT_PRELOAD)

# UI
ui <- dashboardPage(
  skin = "black",

  dashboardHeader(
    title = tags$span(
      style = "display:inline-flex; align-items:center; gap:8px;",
      HTML('<svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="#1a8a6a" stroke-width="2" stroke-linecap="round" style="vertical-align:middle;"><path d="M5 3c0 6 14 6 14 12M5 9c0 6 14 6 14 12"/><path d="M5 3c0 6 14 6 14 12M5 9c0 6 14 6 14 12" transform="scale(-1,1) translate(-24,0)"/><line x1="7" y1="5" x2="17" y2="5"/><line x1="6" y1="9" x2="18" y2="9"/><line x1="6" y1="15" x2="18" y2="15"/><line x1="7" y1="19" x2="17" y2="19"/></svg>'),
      tags$span("HAPView", style = "color:#1a3a4a; font-weight:700;")
    ),
    titleWidth = 220
  ),

  dashboardSidebar(
    width = 220,
    tags$head(
      tags$style(HTML("
        /* ── light theme ── */
        /* accent = teal #1a8a6a ; surfaces = light grey/white ; text = dark slate */
        body, .content-wrapper, .main-footer { background:#eef1f5 !important; }
        .skin-black .main-header .logo { background:#ffffff !important; color:#1a3a4a !important; font-family:'Space Mono',monospace; letter-spacing:1px; border-bottom:1px solid #d7dde6; }
        .skin-black .main-header .logo:hover { background:#f4f6f9 !important; }
        .skin-black .main-header .navbar { background:#ffffff !important; border-bottom:1px solid #d7dde6; }
        .skin-black .main-header .navbar .sidebar-toggle { color:#5a6a7a !important; }
        .skin-black .main-sidebar { background:#f4f6f9 !important; border-right:1px solid #d7dde6; }
        .sidebar-menu > li > a { color:#3a4a5a !important; font-family:'Space Mono',monospace; font-size:12px; }
        .sidebar-menu > li.active > a, .sidebar-menu > li > a:hover { color:#1a8a6a !important; background:#e3eae8 !important; }
        .box { background:#ffffff !important; border:1px solid #d7dde6 !important; border-radius:10px !important; box-shadow:0 1px 3px rgba(0,0,0,0.04); }
        .box-header { background:#f7f9fb !important; border-radius:10px 10px 0 0 !important; color:#1a7a5e !important; border-bottom:1px solid #e3e8ee; }
        .box-title { font-family:'Space Mono',monospace !important; font-size:13px !important; color:#1a7a5e !important; }
        h3.box-title { color:#1a7a5e !important; }
        .content-header h1 { font-family:'Space Mono',monospace; color:#1a3a4a; font-size:20px; }
        .selectize-input { background:#ffffff !important; border:1px solid #c3ccd6 !important; color:#2a3a4a !important; font-size:12px !important; border-radius:6px !important; }
        .selectize-input.focus { border-color:#1a8a6a !important; box-shadow:0 0 0 2px rgba(26,138,106,0.15) !important; }
        .selectize-dropdown { background:#ffffff !important; color:#2a3a4a !important; border:1px solid #c3ccd6 !important; }
        .selectize-dropdown-content .option:hover { background:#e3eae8 !important; }
        .form-control { background:#ffffff !important; border:1px solid #c3ccd6 !important; color:#2a3a4a !important; border-radius:6px !important; font-family:'Space Mono',monospace; font-size:12px; }
        .form-control:focus { border-color:#1a8a6a !important; box-shadow:0 0 0 2px rgba(26,138,106,0.15) !important; }
        .btn-primary { background:#1a8a6a !important; border:none !important; color:#ffffff !important; font-family:'Space Mono',monospace; font-weight:700; border-radius:6px !important; }
        .btn-primary:hover { background:#147055 !important; }
        .btn-default { background:#ffffff !important; border:1px solid #c3ccd6 !important; color:#3a4a5a !important; }
        .shiny-notification { background:#ffffff; color:#1a3a4a; border:1px solid #1a8a6a; font-family:'Space Mono',monospace; box-shadow:0 2px 8px rgba(0,0,0,0.12); }
        label { color:#3a4a5a !important; font-family:'Space Mono',monospace !important; font-size:11px !important; }
        .dataTables_wrapper { color:#3a4a5a; font-family:'Space Mono',monospace; font-size:11px; }
        table.dataTable thead th { color:#1a7a5e !important; background:#f7f9fb !important; border-bottom:1px solid #d7dde6 !important; }
        table.dataTable tbody tr { background:#ffffff !important; color:#2a3a4a !important; }
        table.dataTable tbody tr:hover { background:#f0f4f8 !important; }
        .dataTables_filter input, .dataTables_length select { background:#ffffff !important; border:1px solid #c3ccd6 !important; color:#2a3a4a !important; border-radius:4px; }
        .page-link, .paginate_button { color:#1a8a6a !important; }
        /* upload box */
        .upload-zone { border:2px dashed #c3ccd6; border-radius:12px; padding:30px; text-align:center; color:#7a8a9a; transition:border-color .3s; }
        .upload-zone:hover { border-color:#1a8a6a; }
        /* stat cards */
        .stat-card { background:#ffffff; border:1px solid #d7dde6; border-radius:8px; padding:14px 18px; text-align:center; box-shadow:0 1px 3px rgba(0,0,0,0.04); }
        .stat-num { font-family:'Space Mono',monospace; font-size:22px; font-weight:700; color:#1a8a6a; }
        .stat-lbl { font-family:'Space Mono',monospace; font-size:10px; color:#7a8a9a; text-transform:uppercase; letter-spacing:1px; }
        /* tab styling */
        .nav-tabs { border-bottom:1px solid #d7dde6; }
        .nav-tabs > li > a { color:#5a6a7a !important; font-family:'Space Mono',monospace; font-size:12px; background:#f4f6f9 !important; border:1px solid #d7dde6 !important; }
        .nav-tabs > li.active > a { color:#1a8a6a !important; background:#ffffff !important; border-bottom-color:#ffffff !important; font-weight:700; }
        .nav-tabs > li > a:hover { background:#e9eef3 !important; }
        /* slider */
        .irs--shiny .irs-bar { background:#1a8a6a; }
        .irs--shiny .irs-handle { background:#1a8a6a; border-color:#1a8a6a; }
        .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single { background:#1a8a6a; color:#ffffff; }
        .irs--shiny .irs-line { background:#d7dde6; }
        .irs--shiny .irs-min, .irs--shiny .irs-max { color:#7a8a9a; background:#e3e8ee; }
        /* checkbox labels readable */
        .checkbox label, .radio label { color:#3a4a5a !important; }
        /* file input button */
        .btn-file { background:#1a8a6a !important; color:#fff !important; border:none !important; }
        .progress-bar { background-color:#1a8a6a !important; }
      "))
    ),
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(href = "https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&display=swap", rel = "stylesheet"),

    br(),
    div(style = "padding:0 15px;",
      # Hide file picker entirely when a pre-loaded object is configured
      if (!PRELOAD_MODE) tagList(
        # ── Primary: server-side path (instant, no upload) ──
        textInput("rds_path", "Server file path",
          placeholder = "path/to/object.rds"),
        actionButton("load_path", "Load from path", class = "btn-primary",
          style = "width:100%; margin-bottom:8px;"),
        tags$p(style = "color:#9aa7b4; font-family:'Space Mono',monospace;
                        font-size:10px; text-align:center; margin:6px 0;",
               "── or upload from browser ──"),
        # ── Fallback: browser upload ──
        fileInput("rds_file", NULL,
          accept = c(".rds", ".RDS"),
          buttonLabel = "Browse",
          placeholder = "No file selected"),
        hr(style = "border-color:#d7dde6;")
      ) else tagList(
        tags$p(style = "color:#1a8a6a; font-family:'Space Mono',monospace;
                        font-size:11px; text-align:center; padding:8px 0;",
               "Tasosoid dataset loaded"),
        hr(style = "border-color:#d7dde6;")
      ),
      uiOutput("metadata_col_ui"),
      uiOutput("gene_select_ui"),
      hr(style = "border-color:#d7dde6;"),
      uiOutput("reduction_ui"),
      uiOutput("assay_ui"),
      uiOutput("slot_ui")
    )
  ),

  dashboardBody(
    # Dataset overview stats
    uiOutput("stats_row"),

    tabBox(
      width = 12,
      id = "main_tabs",

      # Tab 1: Dimensionality Reduction 
      tabPanel("UMAP / PCA",
        fluidRow(
          column(6,
            box(width = NULL, title = "Colored by Cell Type",
              plotlyOutput("dimred_celltype", height = "480px")
            )
          ),
          column(6,
            box(width = NULL, title = "Colored by Gene Expression",
              plotlyOutput("dimred_expr", height = "480px")
            )
          )
        ),
        fluidRow(
          column(4,
            box(width = NULL, title = "Plot Options",
              sliderInput("pt_size", "Point size", 0.2, 5, 1.2, step = 0.1),
              sliderInput("pt_alpha", "Opacity", 0.1, 1, 0.8, step = 0.05),
              checkboxInput("show_labels", "Show cluster labels", TRUE),
              checkboxInput("use_custom_colors", "Use atlas color palette", TRUE),
              selectInput("color_scale", "Expression color scale",
                choices = c("FeaturePlot (grey-blue)" = "featureplot",
                            "viridis", "magma", "plasma", "inferno", "turbo"),
                selected = "featureplot")
            )
          ),
          column(8,
            box(width = NULL, title = "Cell Type Composition",
              plotlyOutput("celltype_bar", height = "250px")
            )
          )
        )
      ),

      # Tab 2: Heatmap 
      tabPanel("Heatmap",
        fluidRow(
          column(3,
            box(width = NULL, title = "Heatmap Options",
              uiOutput("heatmap_genes_ui"),
              numericInput("top_n_markers", "Top N markers per cluster", 5, 1, 20),
              tags$p(style="color:#5a6a7a; font-family:\'Space Mono\',monospace; font-size:10px; margin:10px 0 4px;",
                     "FindAllMarkers() thresholds"),
              numericInput("min_pct", "min.pct", 0.05, 0, 1, step = 0.01),
              numericInput("logfc_threshold", "logfc.threshold", 0.25, 0, 5, step = 0.05),
              numericInput("padj_cutoff", "p_val_adj <", 0.05, 0, 1, step = 0.01),
              numericInput("log2fc_cutoff", "avg_log2FC >", 1, 0, 10, step = 0.1),
              actionButton("run_markers", "Find Markers & Plot", class = "btn-primary"),
              br(), br(),
              selectInput("heatmap_scale", "Scale by",
                choices = c("none", "row", "column"), selected = "row"),
              selectInput("heatmap_palette", "Color palette",
                choices = c("RdBu", "PRGn", "RdYlBu", "Spectral", "viridis"),
                selected = "RdBu")
            )
          ),
          column(9,
            box(width = NULL, title = "Expression Heatmap",
              div(style = "overflow-x:auto;",
                plotOutput("heatmap_plot", height = "600px")
              )
            )
          )
        )
      ),

      # Tab 3: Violin / Feature plots 
      tabPanel("Gene Expression",
        fluidRow(
          column(12,
            box(width = NULL, title = "Violin Plot — Gene across Cell Types",
              plotlyOutput("violin_plot", height = "420px")
            )
          )
        ),
        fluidRow(
          column(12,
            box(width = NULL, title = "Dot Plot — Multiple Genes",
              tags$p(style="color:#5a6a7a; font-family:'Space Mono',monospace; font-size:11px; margin-bottom:6px;",
                     "Type to search, click each gene to add it. Selected genes appear as tags. Then click Update."),
              selectizeInput("dotplot_genes", "Genes for dot plot",
                             choices = NULL, multiple = TRUE, selected = NULL,
                             width = "100%",
                             options = list(placeholder = "Type a gene name...")),
              actionButton("run_dotplot", "Update Dot Plot", class = "btn-primary"),
              br(), br(),
              plotOutput("dot_plot", height = "400px")
            )
          )
        )
      ),

      # Tab 4: Cell Type Annotations 
      tabPanel("Annotations",
        fluidRow(
          column(5,
            box(width = NULL, title = "Cell Type Counts",
              DTOutput("celltype_table")
            )
          ),
          column(7,
            box(width = NULL, title = "Cell Distribution (Pie)",
              plotlyOutput("celltype_pie", height = "420px")
            )
          )
        ),
        fluidRow(
          column(12,
            box(width = NULL, title = "Full Cell Metadata",
              DTOutput("metadata_table")
            )
          )
        )
      ),

      # Tab 5: Dataset Info
      tabPanel("Dataset Info",
        box(width = 12, title = "Seurat Object Summary",
          verbatimTextOutput("seurat_summary")
        )
      )
    )
  )
)

# SERVER 
# Increase max upload size to 10 GB
options(shiny.maxRequestSize = 10 * 1024^3)

server <- function(input, output, session) {

  # Track which source last triggered a load
  load_trigger <- reactiveVal(
    # In preload mode, set the trigger immediately so the object loads on startup
    if (PRELOAD_MODE) list(source = "preload", path = PRELOAD_PATH) else NULL
  )

  # Load from server path button (disabled in preload mode)
  observeEvent(input$load_path, {
    p <- trimws(input$rds_path)
    req(nchar(p) > 0)
    load_trigger(list(source = "path", path = p))
  })

  # Load from browser upload (fires automatically when file is chosen)
  observeEvent(input$rds_file, {
    req(input$rds_file)
    load_trigger(list(source = "upload", path = input$rds_file$datapath))
  })

  # Reactive: load Seurat object from whichever source fired
  seurat_obj <- eventReactive(load_trigger(), {
    trigger <- load_trigger(); req(trigger)
    filepath <- trigger$path

    if (!file.exists(filepath)) {
      showNotification(paste("File not found:", filepath), type = "error", duration = 8)
      return(NULL)
    }

    withProgress(message = "Loading Seurat object...", value = 0.2, {
      tryCatch({
        incProgress(0.3, detail = "Reading RDS...")
        obj <- readRDS(filepath)
        if (!inherits(obj, "Seurat")) {
          showNotification("File does not appear to be a Seurat object.", type = "error")
          return(NULL)
        }
        incProgress(0.5, detail = "Done.")
        showNotification(
          paste0("Loaded: ", ncol(obj), " cells x ", nrow(obj), " genes"),
          type = "message", duration = 4
        )
        obj
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error", duration = 8)
        NULL
      })
    })
  }, ignoreNULL = TRUE)

  # Dynamic UI outputs
  output$metadata_col_ui <- renderUI({
    obj <- seurat_obj(); req(obj)
    meta_cols <- colnames(obj@meta.data)
    # Priority 1: exact match for predicted.celltype_updated
    # Priority 2: columns containing "celltype_updated" or "predicted.celltype"
    # Priority 3: any column matching general annotation keywords
    default_val <- if ("predicted.celltype_updated" %in% meta_cols) {
      "predicted.celltype_updated"
    } else {
      specific <- grep("celltype_updated|predicted.celltype", meta_cols,
                       ignore.case = TRUE, value = TRUE)
      if (length(specific) > 0) specific[1] else {
        general <- grep("celltype|cluster|ident|annot|type|label", meta_cols,
                        ignore.case = TRUE, value = TRUE)
        if (length(general) > 0) general[1] else meta_cols[1]
      }
    }
    selectInput("meta_col", "Cell type / annotation column",
                choices = meta_cols, selected = default_val)
  })

  output$gene_select_ui <- renderUI({
    obj <- seurat_obj(); req(obj)
    genes <- rownames(obj)
    selectizeInput("selected_gene", "Gene to visualize",
                   choices = NULL, selected = NULL,
                   options = list(maxOptions = 50, placeholder = "Type gene name..."))
  })

  # Server-side gene search for performance
  observe({
    obj <- seurat_obj(); req(obj)
    updateSelectizeInput(session, "selected_gene",
                         choices = rownames(obj),
                         server = TRUE)
  })

  output$reduction_ui <- renderUI({
    obj <- seurat_obj(); req(obj)
    reds <- names(obj@reductions)
    selectInput("reduction", "Reduction", choices = reds,
                selected = if ("umap" %in% tolower(reds)) reds[grep("umap", reds, ignore.case=TRUE)[1]] else reds[1])
  })

  output$assay_ui <- renderUI({
    obj <- seurat_obj(); req(obj)
    selectInput("assay", "Assay",
                choices = Assays(obj), selected = DefaultAssay(obj))
  })

  output$slot_ui <- renderUI({
    selectInput("expr_slot", "Expression slot",
                choices = c("data", "counts", "scale.data"), selected = "data")
  })

  output$heatmap_genes_ui <- renderUI({
    obj <- seurat_obj(); req(obj)
    selectizeInput("heatmap_custom_genes", "Custom genes (optional)",
                   choices = NULL, multiple = TRUE, selected = NULL,
                   options = list(maxOptions = 50, placeholder = "Leave empty to auto-find markers"))
  })

  observe({
    obj <- seurat_obj(); req(obj)
    updateSelectizeInput(session, "heatmap_custom_genes",
                         choices = rownames(obj), server = TRUE)
  })

  # Populate dot plot gene choices server-side once an object is loaded
  observeEvent(seurat_obj(), {
    obj <- seurat_obj(); req(obj)
    updateSelectizeInput(session, "dotplot_genes",
                         choices = rownames(obj), selected = character(0),
                         server = TRUE)
  })

  # Stats row 
  output$stats_row <- renderUI({
    obj <- seurat_obj()
    if (is.null(obj)) {
      return(fluidRow(
        column(12, div(style = "text-align:center; padding:60px; color:#7a8a9a; font-family:'Space Mono',monospace;",
          tags$p(style="font-size:14px;", "Enter a server path or upload a .rds file to begin"),
          tags$p(style="font-size:11px; color:#a3b0bd;", "Supports Seurat v4 and v5 objects")
        ))
      ))
    }
    n_cells  <- ncol(obj)
    n_genes  <- nrow(obj)
    n_clust  <- if (!is.null(input$meta_col) && input$meta_col %in% colnames(obj@meta.data))
                  length(unique(obj@meta.data[[input$meta_col]])) else "—"
    n_reds   <- length(names(obj@reductions))

    fluidRow(
      column(3, div(class="stat-card",
        div(class="stat-num", format(n_cells, big.mark=",")),
        div(class="stat-lbl", "Cells"))),
      column(3, div(class="stat-card",
        div(class="stat-num", format(n_genes, big.mark=",")),
        div(class="stat-lbl", "Genes"))),
      column(3, div(class="stat-card",
        div(class="stat-num", n_clust),
        div(class="stat-lbl", "Cell types"))),
      column(3, div(class="stat-card",
        div(class="stat-num", n_reds),
        div(class="stat-lbl", "Reductions")))
    )
  })

  # Helper: get reduction coords 
  get_embedding <- reactive({
    obj <- seurat_obj(); req(obj, input$reduction)
    emb <- Embeddings(obj, reduction = input$reduction)
    as.data.frame(emb[, 1:2]) %>%
      setNames(c("Dim1", "Dim2")) %>%
      rownames_to_column("cell")
  })

  # UMAP/PCA — Cell type 
  output$dimred_celltype <- renderPlotly({
    obj <- seurat_obj(); req(obj, input$meta_col, input$reduction)
    emb <- get_embedding()
    meta <- obj@meta.data %>% rownames_to_column("cell")
    df <- left_join(emb, meta[, c("cell", input$meta_col)], by = "cell")
    colnames(df)[4] <- "CellType"

    ax_label <- toupper(input$reduction)

    # Build color scale: use atlas palette when toggled on, fallback to ggplot default
    color_scale_ct <- if (isTRUE(input$use_custom_colors)) {
      scale_color_manual(values = resolve_colors(df$CellType), na.value = "#aaaaaa")
    } else {
      scale_color_discrete()
    }

    p <- ggplot(df, aes(Dim1, Dim2, color = CellType,
                        text = paste0("Cell: ", cell, "<br>Type: ", CellType))) +
      geom_point(size = input$pt_size, alpha = input$pt_alpha) +
      color_scale_ct +
      coord_fixed() +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(color = "#333333", size = 9,
                                   family = "Space Mono"),
        legend.title = element_text(color = "#111111", size = 9,
                                    family = "Space Mono"),
        axis.title = element_text(color = "#666666", size = 9),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      labs(color = "Cell Type",
           x = paste(ax_label, "1"), y = paste(ax_label, "2")) +
      guides(color = "none")

    if (input$show_labels) {
      centroids <- df %>%
        group_by(CellType) %>%
        summarise(Dim1 = median(Dim1), Dim2 = median(Dim2), .groups = "drop")
      p <- p + geom_text(data = centroids,
                         aes(x = Dim1, y = Dim2, label = CellType),
                         color = "#111111", size = 2.5,
                         fontface = "bold", inherit.aes = FALSE)
    }

    ggplotly(p, tooltip = "text") %>%
      layout(
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
        font = list(color = "#333333", family = "Space Mono"),
        showlegend = FALSE
      ) %>%
      config(displayModeBar = FALSE)
  })

  # UMAP/PCA — Gene expression ───────────────────────────
  output$dimred_expr <- renderPlotly({
    obj <- seurat_obj(); req(obj, input$selected_gene, input$reduction)
    gene <- input$selected_gene
    req(gene %in% rownames(obj))

    emb <- get_embedding()
    expr_vals <- FetchData(obj, vars = gene, layer = input$expr_slot)
    df <- cbind(emb, expr = expr_vals[[1]])

    palette_fn <- switch(input$color_scale,
      featureplot = colorRampPalette(c("lightgrey", "blue"))(100),  # Seurat FeaturePlot default
      viridis  = viridis::viridis(100),
      magma    = viridis::magma(100),
      plasma   = viridis::plasma(100),
      inferno  = viridis::inferno(100),
      turbo    = viridis::turbo(100)
    )
    ax_label <- toupper(input$reduction)

    p <- ggplot(df, aes(Dim1, Dim2, color = expr,
                        text = paste0("Cell: ", cell,
                                      "<br>", gene, ": ", round(expr, 3)))) +
      geom_point(size = input$pt_size, alpha = input$pt_alpha) +
      scale_color_gradientn(colors = palette_fn,
                            name = paste0(gene, "\n(expr)")) +
      coord_fixed() +
      theme_void() +
      theme(
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text  = element_text(color = "#333333", size = 9,
                                    family = "Space Mono"),
        legend.title = element_text(color = "#111111", size = 9,
                                    family = "Space Mono")
      ) +
      labs(x = paste(ax_label, "1"), y = paste(ax_label, "2"))

    ggplotly(p, tooltip = "text") %>%
      layout(
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
        font = list(color = "#333333", family = "Space Mono")
      ) %>%
      config(displayModeBar = FALSE)
  })

  # Cell type bar chart
  output$celltype_bar <- renderPlotly({
    obj <- seurat_obj(); req(obj, input$meta_col)
    counts <- as.data.frame(table(obj@meta.data[[input$meta_col]]))
    colnames(counts) <- c("CellType", "Count")
    counts <- counts %>% arrange(desc(Count))

    p <- plot_ly(counts, x = ~reorder(CellType, -Count), y = ~Count,
                 type = "bar",
                 marker = list(color = ~Count,
                               colorscale = "Viridis",
                               line = list(color = "white", width = 0.5)),
                 text = ~Count, textposition = "outside",
                 hoverinfo = "x+y") %>%
      layout(
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
        font = list(color = "#333333", family = "Space Mono", size = 10),
        xaxis = list(title = "", tickfont = list(size = 9), showgrid = FALSE,
                     color = "#444444"),
        yaxis = list(title = "# Cells", gridcolor = "#dddddd", color = "#444444"),
        showlegend = FALSE,
        margin = list(b = 100)
      ) %>%
      config(displayModeBar = FALSE)
    p
  })

  # Markers & Heatmap 
  markers_data <- eventReactive(input$run_markers, {
    obj <- seurat_obj(); req(obj, input$meta_col)
    withProgress(message = "Finding markers (this may take a few minutes)...", value = 0.1, {
      Idents(obj) <- obj@meta.data[[input$meta_col]]
      incProgress(0.1, detail = "Setting identities...")
      result <- tryCatch({
        # Mirrors: FindAllMarkers(obj, only.pos = TRUE, min.pct = X, logfc.threshold = Y)
        #          subset(markers, p_val_adj < Z & avg_log2FC > W)
        DefaultAssay(obj) <- input$assay

        # DEBUG LOGGING 
        message("========== FindAllMarkers debug ==========")
        message("Active assay: ", DefaultAssay(obj))
        message("Assay used:   ", input$assay)
        message("Meta column:  ", input$meta_col)
        message("Idents table:")
        print(table(Idents(obj)))
        layer_names <- tryCatch(Layers(obj[[input$assay]]), error = function(e) character(0))
        message("Layers in assay: ", paste(layer_names, collapse = ", "))
        # Peek at a few values from the data layer
        d <- tryCatch(GetAssayData(obj, assay = input$assay, layer = "data"),
                      error = function(e) NULL)
        if (!is.null(d)) {
          message("Data matrix: ", nrow(d), " genes x ", ncol(d), " cells")
          message("Data range: [", round(min(d), 3), ", ", round(max(d), 3), "]")
          message("Nonzero entries: ", sum(d > 0), " (", round(100*sum(d>0)/length(d), 2), "%)")
        }
        ################################################################

        markers <- FindAllMarkers(obj,
                                  assay           = input$assay,
                                  only.pos        = TRUE,
                                  min.pct         = input$min_pct,
                                  logfc.threshold = input$logfc_threshold,
                                  verbose         = TRUE)
        message("FindAllMarkers returned: ",
                if (is.null(markers)) "NULL" else paste(nrow(markers), "rows"))
        if (!is.null(markers) && nrow(markers) > 0) {
          message("Columns: ", paste(colnames(markers), collapse = ", "))
          message("Head:")
          print(utils::head(markers, 3))
        }
        incProgress(0.6, detail = "Filtering by p_val_adj / log2FC...")

        if (is.null(markers) || nrow(markers) == 0) {
          showNotification("FindAllMarkers returned no markers. Try lowering thresholds.", type = "warning", duration = 8)
          return(NULL)
        }

        # Always enforce a p_val_adj cutoff. If the box is empty/NA, fall back to 0.05.
        padj_thr  <- if (is.null(input$padj_cutoff) || is.na(input$padj_cutoff)) 0.05 else input$padj_cutoff
        log2fc_thr <- if (is.null(input$log2fc_cutoff) || is.na(input$log2fc_cutoff)) 0 else input$log2fc_cutoff

        markers <- subset(markers,
                          p_val_adj  < padj_thr &
                          avg_log2FC > log2fc_thr)

        incProgress(0.8, detail = "Done.")
        if (nrow(markers) == 0) {
          showNotification("No markers passed the p_val_adj / log2FC filter. Try relaxing cutoffs.",
                           type = "warning", duration = 8)
          return(NULL)
        }

        showNotification(paste0("Found ", nrow(markers), " markers across ",
                                length(unique(markers$cluster)), " clusters."),
                         type = "message", duration = 5)
        markers
      }, error = function(e) {
        showNotification(paste("FindAllMarkers error:", e$message), type = "error", duration = 10)
        NULL
      })
      result
    })
  })

  output$heatmap_plot <- renderPlot({
    obj <- seurat_obj(); req(obj, input$meta_col)
    Idents(obj) <- obj@meta.data[[input$meta_col]]

    tryCatch({
      # Choose genes: custom or from markers
      if (!is.null(input$heatmap_custom_genes) && length(input$heatmap_custom_genes) > 0) {
        genes_use <- input$heatmap_custom_genes
      } else {
        req(markers_data())
        markers <- markers_data()
        req(!is.null(markers) && nrow(markers) > 0)
        genes_use <- markers %>%
          group_by(cluster) %>%
          slice_max(order_by = avg_log2FC, n = input$top_n_markers) %>%
          pull(gene) %>%
          unique()
      }

      genes_use <- genes_use[genes_use %in% rownames(obj)]
      validate(need(length(genes_use) > 0, "No valid genes found. Try finding markers first or enter custom genes."))

      # Build expression matrix (averaged per cluster)
      expr_mat <- as.matrix(GetAssayData(obj, assay = input$assay,
                                         layer = "data")[genes_use, , drop = FALSE])
      cell_types <- obj@meta.data[[input$meta_col]]
      avg_mat <- sapply(unique(cell_types), function(ct) {
        cells <- which(cell_types == ct)
        rowMeans(expr_mat[, cells, drop = FALSE])
      })
      colnames(avg_mat) <- unique(cell_types)

      # Color palette
      pal <- switch(input$heatmap_palette,
        RdBu     = rev(brewer.pal(11, "RdBu")),
        PRGn     = brewer.pal(11, "PRGn"),
        RdYlBu   = rev(brewer.pal(11, "RdYlBu")),
        Spectral = rev(brewer.pal(11, "Spectral")),
        viridis  = viridis(50)
      )

      pheatmap(avg_mat,
               scale        = input$heatmap_scale,
               color        = colorRampPalette(pal)(100),
               border_color = NA,
               fontsize_row = 8,
               fontsize_col = 9,
               angle_col    = 45,
               main         = "Mean Expression per Cell Type",
               silent       = FALSE)

    }, error = function(e) {
      plot.new()
      msg <- paste("Heatmap error:", e$message, sep = "\n")
      wrapped <- paste(strwrap(msg, width = 80), collapse = "\n")
      text(0.5, 0.5, wrapped, col = "red", cex = 0.9, family = "mono")
    })
  }, bg = "white")

  # Violin plot
  output$violin_plot <- renderPlotly({
    obj <- seurat_obj(); req(obj, input$selected_gene, input$meta_col)
    gene <- input$selected_gene
    req(gene %in% rownames(obj))

    expr <- FetchData(obj, vars = c(gene, input$meta_col), layer = input$expr_slot)
    colnames(expr) <- c("Expression", "CellType")

    ct_order <- expr %>%
      group_by(CellType) %>%
      summarise(med = median(Expression), .groups = "drop") %>%
      arrange(desc(med)) %>%
      pull(CellType)

    expr$CellType <- factor(expr$CellType, levels = ct_order)

    # Use atlas palette when the UMAP toggle is on; else default qualitative
    violin_colors <- if (isTRUE(input$use_custom_colors)) {
      resolve_colors(levels(expr$CellType))
    } else {
      setNames(scales::hue_pal()(length(levels(expr$CellType))), levels(expr$CellType))
    }

    p <- plot_ly(expr, y = ~Expression, x = ~CellType,
                 type = "violin",
                 box    = list(visible = TRUE),
                 points = "none",
                 color  = ~CellType,
                 colors = violin_colors,
                 hoverinfo = "y+x") %>%
      layout(
        title = list(text = paste0("<b>", gene, "</b>"), font = list(color="#333333")),
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
        font = list(color = "#333333", family = "Space Mono", size = 10),
        xaxis = list(title = "", tickangle = -40, gridcolor = "#eeeeee",
                     color = "#444444"),
        yaxis = list(title = "Expression", gridcolor = "#dddddd", color = "#444444"),
        showlegend = FALSE
      ) %>%
      config(displayModeBar = FALSE)
    p
  })

  # Dot plot
  # Re-render only when "Update Dot Plot" is clicked (genes captured at click time)
  dotplot_genes_active <- eventReactive(input$run_dotplot, {
    input$dotplot_genes
  }, ignoreNULL = FALSE)

  output$dot_plot <- renderPlot({
    obj <- seurat_obj(); req(obj)
    genes <- dotplot_genes_active()
    validate(need(!is.null(genes) && length(genes) > 0,
                  "Select one or more genes above, then click 'Update Dot Plot'."))
    genes <- genes[genes %in% rownames(obj)]
    validate(need(length(genes) > 0, "None of the selected genes were found in this object."))

    Idents(obj) <- obj@meta.data[[input$meta_col]]
    tryCatch({
      DotPlot(obj, features = genes, assay = input$assay) +
        scale_color_gradientn(colors = colorRampPalette(c("lightgrey", "blue"))(100)) +
        coord_flip() +
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.grid       = element_line(color = "#eeeeee"),
          axis.text        = element_text(color = "#333333", size = 9,
                                          family = "Space Mono"),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          axis.title       = element_text(color = "#444444"),
          legend.background = element_rect(fill = "white"),
          legend.text       = element_text(color = "#333333", size = 8,
                                           family = "Space Mono"),
          legend.title      = element_text(color = "#111111", size = 9,
                                           family = "Space Mono")
        )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Dot plot error:\n", e$message), col = "red", cex = 1.0)
    })
  }, bg = "white")

  # Annotations tab
  output$celltype_table <- renderDT({
    obj <- seurat_obj(); req(obj, input$meta_col)
    counts <- as.data.frame(table(obj@meta.data[[input$meta_col]]))
    colnames(counts) <- c("Cell Type", "Count")
    counts$"%" <- paste0(round(counts$Count / sum(counts$Count) * 100, 1), "%")
    counts <- counts %>% arrange(desc(Count))
    datatable(counts, options = list(pageLength = 15, dom = "tp"),
              rownames = FALSE)
  })

  output$celltype_pie <- renderPlotly({
    obj <- seurat_obj(); req(obj, input$meta_col)
    counts <- as.data.frame(table(obj@meta.data[[input$meta_col]]))
    colnames(counts) <- c("CellType", "Count")

    # Match slice colors to the atlas palette when toggle is on
    pie_cols <- if (isTRUE(input$use_custom_colors)) {
      unname(resolve_colors(as.character(counts$CellType))[as.character(counts$CellType)])
    } else {
      scales::hue_pal()(nrow(counts))
    }

    plot_ly(counts, labels = ~CellType, values = ~Count,
            type = "pie",
            marker = list(colors = pie_cols, line = list(color = "white", width = 1.5)),
            textinfo = "label+percent",
            textfont = list(size = 10, family = "Space Mono", color = "#333333"),
            hoverinfo = "label+value") %>%
      layout(
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
        font = list(color = "#333333", family = "Space Mono"),
        legend = list(font = list(color = "#333333", size = 10)),
        margin = list(t = 10, b = 10)
      ) %>%
      config(displayModeBar = FALSE)
  })

  output$metadata_table <- renderDT({
    obj <- seurat_obj(); req(obj)
    meta <- obj@meta.data
    # Show first 500 cells max for performance
    if (nrow(meta) > 500) meta <- meta[1:500, ]
    datatable(meta,
              options = list(pageLength = 10, scrollX = TRUE, dom = "ftip"),
              rownames = TRUE)
  })

  # Dataset summary 
  output$seurat_summary <- renderPrint({
    obj <- seurat_obj(); req(obj)
    print(obj)
    cat("\n── Reductions ──────────────────────────────────\n")
    for (r in names(obj@reductions)) {
      emb <- Embeddings(obj, r)
      cat(sprintf("  %-12s  %d dims\n", r, ncol(emb)))
    }
    cat("\n── Assays ──────────────────────────────────────\n")
    for (a in Assays(obj)) {
      cat(sprintf("  %-12s  %d genes\n", a, nrow(obj[[a]])))
    }
    cat("\n── Meta.data columns ───────────────────────────\n")
    for (col in colnames(obj@meta.data)) {
      vals <- obj@meta.data[[col]]
      info <- if (is.numeric(vals)) sprintf("numeric  [%.2f – %.2f]", min(vals,na.rm=T), max(vals,na.rm=T))
              else sprintf("factor/char  %d unique values", length(unique(vals)))
      cat(sprintf("  %-35s  %s\n", col, info))
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)