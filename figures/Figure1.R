library(ggplot2)
library(dplyr)
library(ggrastr)
#set.seed(42)
source("scripts/TOME_colors.r")

sd_dir  <- "source_data"
out_dir <- "images"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pt_size <- 1.5

tab   <- read.delim(file.path(sd_dir, "TOME_reference_UMAP.tsv.gz"))   # TOME reference
query <- read.delim(file.path(sd_dir, "TOME_query_UMAP.tsv.gz"))       # HAP + TLS cells


## Fig 1E: Qiu 2022 reference UMAP
p <- ggplot(tab, aes(umap_1, umap_2, colour = cell_state)) +
  geom_point(size = pt_size) +
  scale_colour_manual(values = cell_type_updated_colors) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_classic()+
  coord_fixed() 
ggsave(file.path(out_dir, "Fig1E_TOME_reference_UMAP.pdf"), p, width = 10, height = 7)


## Fig 1F UMAPs: HAP and TLS projected onto the TOME reference
query$cell_type_num <- factor(unname(cell_type_ids[query$cell_state]),
                              levels = unname(cell_type_ids))

tab_rep <- tab[rep(seq_len(nrow(tab)), 2), ]

for (cond in c("HAP", "TLS")) {
  umap_subset <- subset(query, condition == cond)
  centroids <- umap_subset %>%
    group_by(cell_type_num, state_id) %>%
    summarise(refumap_1 = median(refumap_1), refumap_2 = median(refumap_2), .groups = "drop")

  p <- ggplot() +
    rasterise(geom_point(data = tab_rep, aes(umap_1, umap_2),
                         size = 3, colour = "grey82", alpha = 0.5), dpi = 300) +
    geom_point(data = umap_subset, aes(refumap_1, refumap_2, colour = cell_type_num), size = 3) +
    scale_colour_manual(name = "", values = cell_type_final_colors) +
    geom_text(data = centroids, aes(refumap_1, refumap_2, label = state_id),
              size = 6, fontface = "bold", colour = "black") +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    coord_fixed() +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, paste0("Fig1F_UMAP_", cond, "_projected_TOME.pdf")), p, width = 10, height = 7)
}

# cell numbers 
# 2438   
# 1379    

## Fig 1F, stacked bar: cell-state composition TLS vs HAP
comp <- query %>%
  count(condition, cell_state, name = "count") %>%
  group_by(condition) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# legend order as in 03 (order of the palette), with the numbered labels from the figure
present <- intersect(names(cell_type_updated_colors), unique(comp$cell_state))
comp$cell_type_num <- factor(unname(cell_type_ids[comp$cell_state]),
                             levels = unname(cell_type_ids[present]))
comp$condition <- factor(comp$condition, levels = c("TLS", "HAP"))

p <- ggplot(comp, aes(condition, percentage, fill = cell_type_num)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_final_colors) +
  labs(x = "", y = "Percentage of cells", fill = "Cell states") +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16, angle = 45, hjust = 1, colour = "black"),
        axis.text.y  = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"),
        legend.text  = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 16, colour = "black"))
ggsave(file.path(out_dir, "Fig1F_TOME_composition.pdf"), p, width = 10, height = 7)


## Fig 1G: Qiu 2024 reference UMAP
source("scripts/OMG_colors.r")


ref   <- read.delim(file.path(sd_dir, "OMG_reference_UMAP.tsv.gz"))
query <- read.delim(file.path(sd_dir, "OMG_query_UMAP.tsv.gz"))

add_num <- function(df) {
  df$cell_type_num <- factor(unname(celltype_order[df$cell_state]),
                             levels = unname(celltype_order))
  df
}
ref   <- add_num(ref)
query <- add_num(query)

theme_umap <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Fig 1G: OMG reference
centroids <- ref %>% group_by(cell_type_num, state_id) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
p <- ggplot() +
  geom_point(data = ref, aes(umap_1, umap_2, colour = cell_type_num), size = 3) +
  scale_colour_manual(name = "", values = cell_type_colored_numbered) +
  geom_text(data = centroids, aes(umap_1, umap_2, label = state_id),
            size = 6, fontface = "bold", colour = "black") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  coord_fixed() + theme_umap
ggsave(file.path(out_dir, "Fig1G_OMG_reference_UMAP.pdf"), p, width = 10, height = 7)

## Fig 1H: HAP projected onto OMG (also ED 9E: run with cond = "HAPX")
ref_rep <- ref[rep(seq_len(nrow(ref)), 5), ]

plot_projection <- function(cond) {
  sub <- subset(query, condition == cond)
  centroids <- sub %>% group_by(cell_type_num, state_id) %>%
    summarise(refumap_1 = median(refumap_1), refumap_2 = median(refumap_2), .groups = "drop")
  ggplot() +
    rasterise(geom_point(data = ref_rep, aes(umap_1, umap_2),
                         size = 3, colour = "grey82", alpha = 0.5), dpi = 300) +
    geom_point(data = sub, aes(refumap_1, refumap_2, colour = cell_type_num), size = 3) +
    scale_colour_manual(name = "", values = cell_type_colored_numbered) +
    geom_text(data = centroids, aes(refumap_1, refumap_2, label = state_id),
              size = 6, fontface = "bold", colour = "black") +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    coord_fixed() + theme_umap
}
ggsave(file.path(out_dir, "Fig1H_UMAP_HAP_projected_OMG.pdf"), plot_projection("HAP"), width = 10, height = 7)

## Fig 1H: composition bar (same function for 4E, 5E, ED 9E with other conditions)
plot_composition <- function(conds, labels = conds) {
  comp <- query %>%
    filter(condition %in% conds) %>%
    count(condition, cell_state, name = "count") %>%
    group_by(condition) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
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
ggsave(file.path(out_dir, "Fig1H_OMG_composition.pdf"), plot_composition("HAP"), width = 7, height = 7)

#1377 cells for HAP 

## Source Data: one sheet per panel
library(openxlsx)
   wb <- createWorkbook()
   addWorksheet(wb, "Fig1E");       writeData(wb, "Fig1E", tab)
   addWorksheet(wb, "Fig1F_UMAP");  writeData(wb, "Fig1F_UMAP", read.delim(file.path(sd_dir, "TOME_query_UMAP.tsv.gz")))
   addWorksheet(wb, "Fig1F_bar");   writeData(wb, "Fig1F_bar", comp)
   addWorksheet(wb, "Fig1G");       writeData(wb, "Fig1G", ref)
   addWorksheet(wb, "Fig1H_UMAP");  writeData(wb, "Fig1H_UMAP", subset(query, condition == "HAP"))
   addWorksheet(wb, "Fig1H_bar");   writeData(wb, "Fig1H_bar", subset(query, condition == "HAP") |>
                                                 dplyr::count(cell_state) |>
                                                 dplyr::mutate(percentage = n / sum(n) * 100))
   saveWorkbook(wb, file.path(sd_dir, "Source_Data_Fig1.xlsx"), overwrite = TRUE)