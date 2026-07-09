##############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

cell_type_colors = c(
## Neural ectoderm  (brown-orange)
"Anterior floor plate"  =  "#a14a57", 
"Posterior floor plate"  = "#a0522b",          
"Hindbrain" = "#ff8c00",            
"Mesencephalon/MHB" = "#ffae42",
"Neuromesodermal progenitors" = "#9b1c31",
"Spinal cord (dorsal)" = "#e3b778",      
"Spinal cord" = "#ffdcad", 
"Motor neurons" = "#dc143c",
"Rostral neuroectoderm" = "#ff7f50",
"Di/telencephalon" = "#ff6447",               
"Forebrain/midbrain" = "#c1996b",


## Neural 
"Neuron progenitor cells" = "#E38C8C",    
"Spinal cord excitatory neurons" = "#DB6B6B",

## Mesoderm 
"Endothelium" = "#d1f1d2",
"Intermediate mesoderm" = "#4f4a83",    
"Paraxial mesoderm A" = "#97e2ff",       
"Paraxial mesoderm B" = "#156888",
"Somatic mesoderm" = "#1d3557",


# Endo
"Gut" =  "#f6f201",
"Definitive endoderm" = "#ffd700",             

## Intestine 
"Extraembryonic visceral endoderm" =  "#f6c100",

"Skeletal muscle progenitors"  = "#cd5c5c",  
"Hematoendothelial progenitors" = "#6c8e23",   
"White blood cells" = "#228B22",


## epithelial
 "Placodal area" = "#DDA0DD",

"Primordial germ cells" =  "#ff8dfc",       

"Splanchnic mesoderm" = "#c600cb",


"Extraembryonic ectoderm" = "#ffa07a",         
"Pre-epidermal keratinocytes" = "#ffb6c1")


day_colors_f = c("E7.25"="#7b201f", "E7.5"="#e07a5f", "E7.75" = "#f2cc8f", "E8"="#e2953c", "E8.25"="#80b29a", "E8.5"= "#78c1e2",  "E9.5"="#3e405b", "E10.5"="#673e91")
day_colors = c("E7.25"="#7b201f", "E7.5"="#e07a5f", "E7.75" = "#f2cc8f", "E8.0"="#e2953c", "E8.25"="#80b29a", "E8.5"= "#78c1e2",  "E9.5"="#3e405b", "E10.5"="#673e91")


cell_type_updated_colors = c(
## Neural ectoderm  (brown-orange)
"Anterior floor plate"  =  "#a14a57", 
"Posterior floor plate"  = "#a0522b",          
"Hindbrain" = "#ff8c00",  
"Mesencephalon/MHB" = "#ffae42",
"Neuromesodermal progenitors" = "#9b1c31",
"Spinal cord (dorsal)" = "#e3b778",      
"Spinal cord" = "#ffdcad", 
"Motor neurons" = "#dc143c",
"Rostral neuroectoderm" = "#ff7f50",
"Di/telencephalon" = "#ff6447",               
"Midbrain" = "#c1996b",
"Forebrain" = "#ff4500",
"Unassigned 1" = "#927554",


## Neural 
 "Neuron progenitor cells" = "#E38C8C",
 "Unassigned 2" = "#783B3B",
"Spinal cord excitatory neurons" = "#DB6B6B",

## Mesoderm 
"Endothelium" = "#d1f1d2",
"Intermediate mesoderm" = "#4f4a83",    
"Paraxial mesoderm A" = "#97e2ff",       
"Paraxial mesoderm B" = "#156888",
"Somatic mesoderm" = "#1d3557",


# Endo
"Gut" =  "#f6f201",
"Definitive endoderm" = "#ffd700",             

## Intestine 
"Extraembryonic visceral endoderm" =  "#f6c100",

"Skeletal muscle progenitors"  = "#8a8fff",  
"Hematoendothelial progenitors" = "#6c8e23",   
"White blood cells" = "#228B22",


## epithelial
 "Placodal area" = "#DDA0DD",

"Primordial germ cells" =  "#ff8dfc",       

"Splanchnic mesoderm" = "#c600cb",


"Extraembryonic ectoderm" = "#ffa07a",         
"Pre-epidermal keratinocytes" = "#ffb6c1"
)


cell_type_ids <- c(
  "Neuromesodermal progenitors"   = "Neuromesodermal progenitors (1)",
  "Forebrain"                     = "Forebrain (2)",
  "Midbrain"                      = "Midbrain (3)",
  "Hindbrain"                     = "Hindbrain (4)",
  "Spinal cord"                   = "Spinal cord (5)",
  "Paraxial mesoderm A"           = "Paraxial mesoderm A (6)",
  "Paraxial mesoderm B"           = "Paraxial mesoderm B (7)",
  "Hematoendothelial progenitors" = "Hematoendothelial progenitors (8)",
  "Endothelium"                   = "Endothelium (9)",
  "Gut"                           = "Gut (10)"
)

cell_type_final_colors = c(
  "Neuromesodermal progenitors (1)"   = "#9b1c31",
  "Forebrain (2)"                     = "#ff4500",
  "Midbrain (3)"                      = "#c1996b",
  "Hindbrain (4)"                     = "#ff8c00",
  "Spinal cord (5)"                   = "#ffdcad",
  "Paraxial mesoderm A (6)"          = "#97e2ff",
  "Paraxial mesoderm B (7)"          = "#156888",
  "Hematoendothelial progenitors (8)"= "#6c8e23",
  "Endothelium (9)"                  = "#d1f1d2",
  "Gut (10)"                          = "#f6f201"
)
