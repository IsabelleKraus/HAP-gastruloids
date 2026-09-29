library(ggplot2)
library(dplyr)
library(openxlsx)
source("scripts/OMG_colors.r")

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hap   <- read.delim(file.path(sd_dir, "HAP_UMAP_expression.tsv.gz"), check.names = FALSE)
nap   <- read.delim(file.path(sd_dir, "NAP_on_HAP_UMAP_expression.tsv.gz"), check.names = FALSE)
query <- read.delim(file.path(sd_dir, "OMG_query_UMAP.tsv.gz"))    # all conditions, for the bar

# numbered labels
add_num <- function(df) {
  df$cell_type_num <- factor(unname(celltype_order[df$cell_state]),
                             levels = unname(celltype_order))
  df
}
hap <- add_num(hap)
nap <- add_num(nap)

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# UMAP coloured and numbered by OMG cell state; bg = optional grey background
# (06: single grey HAP layer, alpha 0.5, not rasterised)
plot_states <- function(df, bg = NULL) {
  centroids <- df %>% filter(!is.na(cell_type_num)) %>%
    group_by(cell_type_num, state_id) %>%
    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
  p <- ggplot()
  if (!is.null(bg)) p <- p + geom_point(data = bg, aes(umap_1, umap_2),
                                        colour = "grey82", size = 3, alpha = 0.5)
  p +
    geom_point(data = df, aes(umap_1, umap_2, colour = cell_type_num), size = 3) +
    scale_colour_manual(name = "", values = cell_type_colored_numbered, na.value = "gray80") +
    geom_text(data = centroids, aes(umap_1, umap_2, label = state_id),
              size = 6, fontface = "bold", colour = "black") +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    coord_fixed() + theme_umap
}

## Fig 4E: HAP UMAP  
ggsave(file.path(out_dir, "Fig4E_HAP_UMAP.pdf"), plot_states(hap), width = 10, height = 7)

## Fig 4E: NAP projected onto the HAP UMAP (06: UMAP_Normoxic_on_HAP.pdf)
ggsave(file.path(out_dir, "Fig4E_NAP_on_HAP_UMAP.pdf"), plot_states(nap, bg = hap), width = 10, height = 7)

## Fig 4E: composition bar HAP vs NAP (04: stacked_OMG_predicted_states_conditions.pdf)
##   percentage of cells per predicted OMG cell state within each condition
comp_table <- function(conds) {
  query %>%
    filter(condition %in% conds) %>%
    count(condition, cell_state, name = "count") %>%
    group_by(condition) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
}
plot_composition <- function(conds, labels = conds) {
  comp <- comp_table(conds)
  present <- intersect(names(cell_type_colored), unique(comp$cell_state))   # legend order as in 04
  comp$cell_type_num <- factor(unname(celltype_order[comp$cell_state]),
                               levels = unname(celltype_order[present]))
  comp$condition <- factor(comp$condition, levels = conds, labels = labels)
  ggplot(comp, aes(condition, percentage, fill = cell_type_num)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cell_type_colored_numbered) +
    labs(x = "", y = "Percentage of cells", fill = "Cell states") +
    theme_classic() +
    theme(axis.text.x  = element_text(size = 16, angle = 45, hjust = 1, colour = "black"),
          axis.text.y  = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black"),
          legend.text  = element_text(size = 14, colour = "black"),
          legend.title = element_text(size = 16, colour = "black"))
}
ggsave(file.path(out_dir, "Fig4E_OMG_composition.pdf"), plot_composition(c("HAP", "NAP")),
       width = 9, height = 7)


## Source Data
xlsx_file <- file.path(sd_dir, "Source_data_Figure_4.xlsx")
if (file.exists(xlsx_file)) file.copy(xlsx_file, sub("\\.xlsx$", "_backup.xlsx", xlsx_file), overwrite = TRUE)
wb <- if (file.exists(xlsx_file)) loadWorkbook(xlsx_file) else createWorkbook()

put_sheet <- function(wb, name, df) {         # (re)write one sheet
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, df)
}
put_sheet(wb, "Figure 4E HAP UMAP",   hap[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
put_sheet(wb, "Figure 4E NAP on HAP", nap[, c("cell_id", "umap_1", "umap_2", "cell_state", "state_id")])
put_sheet(wb, "Figure 4E bar",        comp_table(c("HAP", "NAP")))
saveWorkbook(wb, xlsx_file, overwrite = TRUE)

