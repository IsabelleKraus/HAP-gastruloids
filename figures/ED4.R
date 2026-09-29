library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
source("scripts/OMG_colors.r")    # cell_type_colored, celltype_order
source("scripts/TOME_colors.r")   # cell_type_ids, cell_type_updated_colors

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Source Data workbook: opened once here, saved once at the end
xlsx_file <- file.path(sd_dir, "Source_data_Extended_Data_Figure_4.xlsx")
if (file.exists(xlsx_file)) file.copy(xlsx_file, sub("\\.xlsx$", "_backup.xlsx", xlsx_file), overwrite = TRUE)
wb <- if (file.exists(xlsx_file)) loadWorkbook(xlsx_file) else createWorkbook()

put_sheet <- function(wb, name, df) {         # (re)write one sheet
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df)
}


## ED Fig 4A: fraction of HAP cells per predicted mouse stage
##   dot size = % of HAP cells, no dot = 0 cells

theme_dots <- function(aspect) {
  theme_classic() +
    theme(axis.text.x  = element_text(size = 14, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y  = element_text(size = 14, colour = "black"),
          axis.title   = element_text(size = 16),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 14),
          aspect.ratio = aspect)
}

## left: Qiu et al. 2022 (TOME); E8.5a / E8.5b of the reference fused into E8.5
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
stage_tome_hap <- filter(stage_tome, condition == "HAP")          # only HAP is shown

p <- ggplot(stage_tome_hap, aes(x = predicted_stage, y = condition, size = percentage)) +
  geom_point(alpha = 0.8, na.rm = TRUE) +
  scale_size_continuous(range = c(2, 10)) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Predicted stage (E)", y = "", size = "Cell fraction (%)", title = "") +
  theme_dots(3/5)
ggsave(file.path(out_dir, "ED4A_TOME_stage_HAP.pdf"), p, width = 7, height = 3)
put_sheet(wb, "ED Fig 4A Qiu 2022", stage_tome_hap)

## right: Qiu et al. 2024 (OMG)
query <- read.delim(file.path(sd_dir, "OMG_query_UMAP.tsv.gz"))
all_stages <- c("E8.0-E8.5", "E8.75", "E9.0", "E9.25", "E9.5", "E9.75")

stage_omg <- query %>%
  filter(condition == "HAP") %>%
  count(condition, predicted_day, name = "cell_count") %>%
  mutate(percentage = cell_count / sum(cell_count) * 100) %>%
  complete(condition = "HAP", predicted_day = all_stages,
           fill = list(cell_count = NA, percentage = NA))
stage_omg$predicted_day <- factor(stage_omg$predicted_day, levels = all_stages, ordered = TRUE)

p <- ggplot(stage_omg, aes(x = predicted_day, y = condition, size = percentage)) +
  geom_point(alpha = 0.8, na.rm = TRUE) +
  scale_size_continuous(range = c(2, 10)) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Predicted stage (E)", y = "", size = "Cell fraction (%)", title = "") +
  theme_dots(0.5)
ggsave(file.path(out_dir, "ED4A_OMG_stage_HAP.pdf"), p, width = 7, height = 2.8)
put_sheet(wb, "ED Fig 4A Qiu 2024", stage_omg)


## ED Fig 4B: score-weighted predicted stage per HAP cell, by predicted cell state

es <- read.delim(file.path(sd_dir, "HAP_expected_stage_mouse.tsv.gz"))

theme_box <- theme_classic() +
  theme(axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, colour = "black"),
        axis.text.y  = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16))

## top: Qiu et al. 2022 (TOME)
##   note: limits = c(8, 9.5) drops cells outside this range before the boxes are computed
d_tome <- filter(es, reference == "Qiu 2022 (TOME)")
d_tome$cell_state <- factor(d_tome$cell_state, levels = names(cell_type_ids))
p <- ggplot(d_tome, aes(cell_state, expected_stage, fill = cell_state)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(y = "Predicted stage", x = "") +
  scale_fill_manual(values = cell_type_updated_colors) +
  scale_y_continuous(breaks = c(8.0, 8.5, 9.0, 9.5, 10.0), limits = c(8, 9.5)) +
  theme_box + theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 60))
ggsave(file.path(out_dir, "ED4B_TOME_stage_per_state.pdf"), p, width = 12, height = 8)
put_sheet(wb, "ED Fig 4B Qiu 2022", d_tome)

## bottom: Qiu et al. 2024 (OMG)
d_omg <- filter(es, reference == "Qiu 2024 (OMG)")
d_omg$cell_state <- factor(d_omg$cell_state, levels = names(cell_type_colored))
p <- ggplot(d_omg, aes(cell_state, expected_stage, fill = cell_state)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(y = "Predicted stage", x = "") +
  scale_fill_manual(values = cell_type_colored) +
  scale_y_continuous(breaks = c(8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75)) +  # stage midpoints + 8.5, as in 04
  theme_box
ggsave(file.path(out_dir, "ED4B_OMG_stage_per_state.pdf"), p, width = 12, height = 6)
put_sheet(wb, "ED Fig 4B Qiu 2024", d_omg)



## ED Fig 4E: HAP trajectory (Monocle3) on the HAP UMAP

tr <- read.delim(file.path(sd_dir, "HAP_trajectory_cells.tsv.gz"))   # one row per HAP cell
gr <- read.delim(file.path(sd_dir, "HAP_trajectory_graph.tsv.gz"))   # principal-graph edges

# numbered labels, e.g. "Midbrain (3)"; factor levels fix the legend order
tr$cell_type_num <- factor(unname(celltype_order[tr$cell_state]),
                           levels = unname(celltype_order))

# as in the figure: axes with ticks, "UMAP 1" / "UMAP 2"
theme_traj <- theme_classic() +
  theme(axis.text = element_text(size = 10, colour = "black"))

# principal graph as drawn by plot_cells(): segments between graph nodes
graph_layer <- function(colour) {
  geom_segment(data = gr, aes(x = x, y = y, xend = xend, yend = yend),
               colour = colour, linewidth = 0.75)
}

## left: cells coloured by OMG cell state
cent <- tr %>% filter(!is.na(cell_type_num)) %>%
  group_by(cell_type_num, state_id) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
p <- ggplot(tr, aes(umap_1, umap_2)) +
  geom_point(aes(colour = cell_type_num), size = 3) +
  graph_layer("black") +
  geom_text(data = cent, aes(label = state_id), size = 5, fontface = "bold") +
  scale_colour_manual(values = cell_type_colored_numbered, name = "") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "") +
  coord_fixed() + theme_traj
ggsave(file.path(out_dir, "ED4E_HAP_trajectory_states.pdf"), p, width = 10, height = 7)

## middle: cells split by the fork-node neighbour their path leaves through;
##   early cells drawn on top
tr_early_on_top <- tr[order(ifelse(is.na(tr$pseudotime), Inf, tr$pseudotime), decreasing = TRUE), ]
p <- ggplot(tr_early_on_top, aes(umap_1, umap_2)) +
  geom_point(aes(colour = root_direction), size = 3) +
  graph_layer("grey28") +                         # plot_cells() default graph colour
  scale_colour_manual(name   = "root_direction",
                      values = c(direction_Brain = "#EAA448", direction_Somites = "#92B9BD",
                                 root = "#9B1C31"),
                      breaks = c("direction_Brain", "direction_Somites", "root"),
                      labels = c("direction_neural", "direction_mesodermal", "root (NMPs)")) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "") +
  coord_fixed() + theme_traj
ggsave(file.path(out_dir, "ED4E_HAP_trajectory_root_directions.pdf"), p, width = 8, height = 6)

## right: Monocle3 pseudotime; grey = cells not connected to the trajectory
p <- ggplot(tr, aes(umap_1, umap_2)) +
  geom_point(aes(colour = pseudotime), size = 3) +
  graph_layer("black") +
  scale_colour_viridis_c(option = "plasma", name = "pseudotime", na.value = "grey80") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "") +
  coord_fixed() + theme_traj
ggsave(file.path(out_dir, "ED4E_HAP_trajectory_pseudotime.pdf"), p, width = 8, height = 6)

put_sheet(wb, "ED Fig 4E cells", tr[, c("cell_id", "umap_1", "umap_2", "cell_state",
                                        "state_id", "pseudotime", "root_direction")])
put_sheet(wb, "ED Fig 4E trajectory graph", gr)

saveWorkbook(wb, xlsx_file, overwrite = TRUE)