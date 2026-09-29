library(ggplot2)
library(dplyr)
library(openxlsx)
source("scripts/OMG_colors.r")

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# TRUE : colour scale 0 - max(HAP, HAP + XAV) per gene 
shared_scale <- FALSE

hapx <- read.delim(file.path(sd_dir, "HAPX_UMAP_expression.tsv.gz"), check.names = FALSE)
vmax <- read.delim(file.path(sd_dir, "HAP_HAPX_shared_vmax.tsv.gz"))
vmax <- setNames(vmax$vmax, vmax$gene)

# numbered labels, e.g. "Anterior Forebrain (1)"; factor levels fix the legend order
hapx$cell_type_num <- factor(unname(celltype_order[hapx$cell_state]),
                             levels = unname(celltype_order))

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Fig 7F: HAP + XAV UMAP, coloured and numbered by OMG cell state
centroids <- hapx %>% filter(!is.na(cell_type_num)) %>%
  group_by(cell_type_num, state_id) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
p <- ggplot(hapx, aes(umap_1, umap_2, colour = cell_type_num)) +
  geom_point(size = 3) +
  scale_colour_manual(name = "", values = cell_type_colored_numbered, na.value = "gray80") +
  geom_text(data = centroids, aes(umap_1, umap_2, label = state_id),
            size = 6, fontface = "bold", colour = "black") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  coord_fixed() + theme_umap
ggsave(file.path(out_dir, "Fig7F_HAPX_UMAP.pdf"), p, width = 10, height = 7)


## Fig 7G: marker expression on the HAP + XAV UMAP
plot_feature <- function(df, gene, limit = NULL) {
  scale <- if (is.null(limit)) {
    scale_colour_gradientn(colours = c("lightgrey", "blue"), name = NULL)
  } else {
    scale_colour_gradientn(colours = c("lightgrey", "blue"), name = NULL,
                           limits = c(0, limit), oob = scales::squish)
  }
  ggplot(df, aes(umap_1, umap_2, colour = .data[[gene]])) +
    geom_point(size = 3) +
    scale +
    ggtitle(gene) + coord_fixed() + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
}
genes_7G <- c("Emx2", "En1", "Lhx2", "Dmbx1",            # brain
              "Nkx6-1", "Olig2", "Pax3", "Pax7",         # neural tube
              "Sox9", "Pax9",                            # sclerotome
              "Foxa2", "Shh",                            # floor plate / notochord
              "Tbx18", "Uncx")                           # somitic (Uncx = Uncx4.1)
for (g in genes_7G) {
  lim <- if (shared_scale) vmax[[g]] else NULL
  ggsave(file.path(out_dir, paste0("Fig7G_", g, "_HAPX_FeaturePlot.pdf")),
         plot_feature(hapx, g, lim), width = 10, height = 7)
}


## Source Data
xlsx_file <- file.path(sd_dir, "Source_data_Figure_7.xlsx")
if (file.exists(xlsx_file)) file.copy(xlsx_file, sub("\\.xlsx$", "_backup.xlsx", xlsx_file), overwrite = TRUE)
wb <- if (file.exists(xlsx_file)) loadWorkbook(xlsx_file) else createWorkbook()

put_sheet <- function(wb, name, df) {         # (re)write one sheet
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df)
}
put_sheet(wb, "Figure 7F", hapx[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
put_sheet(wb, "Figure 7G", hapx[, c("cell_id", "umap_1", "umap_2", "cell_state", genes_7G)])
if (shared_scale) {
  put_sheet(wb, "Figure 7G colour scale",
            data.frame(gene = genes_7G, scale_min = 0, scale_max = unname(vmax[genes_7G])))
}
saveWorkbook(wb, xlsx_file, overwrite = TRUE)

tome_q <- read.delim(file.path(sd_dir, "TOME_query_UMAP.tsv.gz"))


## ED Fig 4A (left): fraction of HAP cells per TOME-predicted stage (Qiu et al. 2022)
##   E8.5a / E8.5b of the reference are fused into E8.5; dot size = % of cells
##   within each condition, no dot = 0 cells
tome_q <- read.delim(file.path(sd_dir, "TOME_query_UMAP.tsv.gz"))
all_stages_tome <- c("E7.25", "E7.5", "E7.75", "E8", "E8.25", "E8.5", "E9.5", "E10.5")

stage_tome <- tome_q %>%
  mutate(predicted_stage = ifelse(predicted_day %in% c("E8.5a", "E8.5b"), "E8.5",
                                  as.character(predicted_day))) %>%
  count(condition, predicted_stage, name = "cell_count") %>%
  group_by(condition) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100) %>%
  ungroup() %>%
  complete(condition = c("TLS", "HAP"), predicted_stage = all_stages_tome,
           fill = list(cell_count = NA, percentage = NA))
stage_tome$predicted_stage <- factor(stage_tome$predicted_stage, levels = all_stages_tome, ordered = TRUE)
stage_tome$condition       <- factor(stage_tome$condition, levels = c("TLS", "HAP"))

p <- ggplot(filter(stage_tome, condition == "HAP"),       # only HAP is shown in ED 4A
            aes(x = predicted_stage, y = condition, size = percentage)) +
  geom_point(alpha = 0.8, na.rm = TRUE) +
  scale_size_continuous(range = c(2, 10)) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Predicted stage (E)", y = "", size = "Cell fraction (%)", title = "") +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(size = 14, colour = "black"),
        axis.title   = element_text(size = 16),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        aspect.ratio = 3/5)
ggsave(file.path(out_dir, "ED4A_TOME_stage_HAP.pdf"), p, width = 7, height = 3)

put_sheet(wb, "ED Fig 4A Qiu 2022", filter(stage_tome, condition == "HAP"))

 