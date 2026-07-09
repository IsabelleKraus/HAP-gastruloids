##############################################################################
# All file paths in this script are relative to the repository root.
# Set the working directory to the repo root before running, e.g.:
#   setwd("/path/to/this/repo")  
# Expected layout: data/  images/  tables/  scripts/  (create as needed)
###############################################################################

cell_type_colors = c(
                    ## Neural ectoderm  (brown-orange)
                    "Anterior floor plate" = "#a14a57",
                    "Floorplate and p3 domain" = "#a0522b",
                    "Anterior roof plate" = "#FFA69E",
                    "Posterior roof plate" = "#ffb6a0",
                    "Spinal cord ventral progenitors" = "#e97451",
                    "Spinal cord motor neurons" = "#ffaf87",
                    "Spinal cord dorsal progenitors" = "#e3b778",
                    "Spinal cord/r7/r8" = "#ffdcad",
                    "NMPs and spinal cord progenitors" = "#9b1c31",
                    "Dorsal telencephalon" = "#ff7f50",
                    "Telencephalon" = "#c50061",
                    "Diencephalon" = "red",
                    "Midbrain" = "#c1996b",
                    "Hindbrain" = "#d94e1f",
                    "Midbrain-hindbrain boundary" = "#ffae42",
                    "Hypothalamus" = "#CC3131",
                    "Hypothalamus (Sim1+)" = "#D34A4A",
                    "Eye field" = "#ff4500",
                    "Cranial motor neurons" = "#dc143c",
                    "Neural crest (PNS neurons)" = "#732507",
                    "Neural crest (PNS glia)" = "#8b4513",

                    ## Surface ectoderm
                    "Olfactory epithelial cells" = "#ff1493",
                    "Otic epithelial cells" = "#ff69b4",
                    "Pre-epidermal keratinocytes" = "#ffb6c1",

                    # Mesoderm
                    # Paraxial mesoderm
                    "Dermomyotome" = "#97e2ff",
                    "Sclerotome" = "#156888",
                    "Mesodermal progenitors (Tbx6+)" = "#6689a1",
                    # Intermediate mesoderm
                    "Anterior intermediate mesoderm" = "#4f4a83",
                    "Lateral plate and intermediate mesoderm" = "#1d3557",
                    # Lateral plate / cardiovascular
                    "First heart field" = "#c600cb",
                    "Second heart field" = "#800080",
                    "Endocardial cells" = "#E3B5E3",
                    "Arterial endothelial cells" = "#da70d6",
                    "Endothelium" = "#d1f1d2",
                    "Pericytes" = "#90ee90",
                    "Hematoendothelial progenitors" = "#6c8e23",
                    "Primitive erythroid cells" = "#556b2f",

                    # Other mesodermal
                    "Facial mesenchyme" = "#4b0082", 
                    "Chondrocytes (Atp1a2+)" = "#8a2be2",
                     "Limb mesenchyme progenitors" = "#9370db",


                    # Endoderm
                    "Gut" = "#f6f201",
                    "Ciliated nodal cells" = "#ffd700",
                    

                    # Germline
                    "Primordial germ cells" = "#ff8dfc",

                    # epithelial
                    "Notochord" = "#1242fd",
                    "Olfactory pit cells" = "#1e90ff",
                    "Granular keratinocytes" = "#87cefa",
                    "Pancreatic acinar cells" = "#00bfff",
                    "Pituitary/Pineal gland progenitors" = "#4682b4",
                    "Apical ectodermal ridge" = "#5f9ea0",
                    "Amniotic ectoderm" = "#ffa07a",
                    "Placodal area" = "#DDA0DD",
                    "Posterior intermediate mesoderm" = "#7B68EE",

                    # Neural crest PNS neurons
                    "Otic sensory neurons" = "#FF69B4",
                    "Sympathetic neurons" = "#FF47A9",
                    "Enteric neurons" = "#A30057",

                    "Olfactory sensory neurons" = "#FF26FF",
               
                    # muscle 
                    "Muscle progenitor cells" = "#cd5c5c",
                    "Muscle progenitor cells (Prdm1+)" = "#b22222",
                    "Myotubes" = "#8b0000",
                    # Hepatocytes
                    "Hepatocytes" = "#daa520",
                    # Lung 
                    "Lung progenitor cells" = "#ffdead",
                    # Eye and other 
                    "Pancreatic islets" = "#20B2AA",
                    "Suprachiasmatic nucleus" = "#008080",
                    # Megakaryocytes 
                    "Megakaryocytes" = "#2E8B57",
                    # Blood 
                    "Hematopoietic stem cells (Cd34+)" = "#006400",
                    "Definitive early erythroblasts (CD36-)" = "#228B22",
                    "Border-associated macrophages" = "#32CD32",
                    "Border-associated macrophages (Ms4a8a+)" = "#3CB371",

                    # CNS neurons
                    "GABAergic neurons" = "#783B3B",
                    "GABAergic cortical interneurons" = "#5C2A2A",
                    "Glutamatergic neurons" = "#DB6B6B",
                    "Neural progenitor cells (Neurod1+)" = "#E38C8C",
                    # Ependymal cells
                    "Choroid plexus" = "#6C350F",
                    # Endothelium 
                    "Brain capillary endothelial cells" = "#AF5931",
                    "Liver sinusoidal endothelial cells" = "#FFA500",
                    
                   # Intestine
                    "Extraembryonic visceral endoderm" = "#f6c100",
                    "Midgut/Hindgut epithelial cells" = "#ffecb3"
)

cell_type_colored = c(
                    ## Neural ectoderm  (brown-orange)
                    "Anterior floor plate" =  "#ff8ba7",
                    "Floorplate and p3 domain" = "#a0522b",
                    "Anterior roof plate" = "#FFA69E",
                    "Posterior roof plate" = "#ffb6a0",
                    "Spinal cord ventral progenitors" = "#e97451",
                    "Spinal cord motor neurons" = "#ffaf87",
                    "Spinal cord dorsal progenitors" = "#e3b778",
                    "Spinal cord/r7/r8" = "#ffdcad",
                    "NMPs and spinal cord progenitors" = "#9b1c31",
                    "Dorsal telencephalon" = "#ff7f50",
                    "Telencephalon" = "#c50061",
                    "Posterior Forebrain / Diencephalon" = "red",
                    "Midbrain" = "#c1996b",
                    "Hindbrain" = "#d94e1f",
                    "Midbrain-hindbrain boundary" = "#ffae42",
                    "Hypothalamus" = "#CC3131",
                    "Hypothalamus (Sim1+)" = "#D34A4A",
                    "Eye field" = "#c50062",
                    "Cranial motor neurons" = "#dc143c",
                    "Neural crest (PNS neurons)" = "#732507",
                    "Neural crest (PNS glia)" = "#8b4513",

                    ## Surface ectoderm
                    "Olfactory epithelial cells" = "#ff1493",
                    "Otic epithelial cells" = "#ff69b4",
                    "Pre-epidermal keratinocytes" = "#ffb6c1",

                    # Mesoderm
                    # Paraxial mesoderm
                    "Somites / Dermomyotome" = "#97e2ff",
                    "Somites / Sclerotome" = "#156888",
                    "Mesodermal progenitors (Tbx6+)" = "#6689a1",
                    # Intermediate mesoderm
                    "Anterior intermediate mesoderm" = "#4f4a83",
                    "Lateral plate and intermediate mesoderm" = "#1d3557",
                    # Lateral plate / cardiovascular
                    "First heart field" = "#c600cb",
                    "Second heart field" = "#800080",
                    "Endocardial cells" = "#E3B5E3",
                    "Arterial endothelial cells" = "#da70d6",
                    "Endothelium" = "#d1f1d2",
                    "Pericytes" = "#90ee90",
                    "Hematoendothelial progenitors" = "#6c8e23",
                    "Primitive erythroid cells" = "#556b2f",

                    # Other mesodermal
                    "FM_0" = "#4b0082", 
                    "FM_1" = "#6e5586",
                    "FM_2" = "#b084c0", 
                    "Chondrocytes (Atp1a2+)" = "#8a2be2",
                     "Limb mesenchyme progenitors" = "#9370db",


                    # Endoderm
                    "Gut" = "#f6f201",
                    "Ciliated nodal cells" = "#ffd700",
                    

                    # Germline
                    "Primordial germ cells" = "#ff8dfc",

                    # epithelial
                    "Notochord" = "#1242fd",
                    "Olfactory pit cells" = "#1e90ff",
                    "Granular keratinocytes" = "#87cefa",
                    "Pancreatic acinar cells" = "#00bfff",
                    "Pituitary/Pineal gland progenitors" = "#4682b4",
                    "Apical ectodermal ridge" = "#5f9ea0",
                    "Amniotic ectoderm" = "#ffa07a",
                    "Placodal area" = "#DDA0DD",
                    "Posterior intermediate mesoderm" = "#7B68EE",

                    # Neural crest PNS neurons
                    "Otic sensory neurons" = "#FF69B4",
                    "Sympathetic neurons" = "#FF47A9",
                    "Enteric neurons" = "#A30057",

                    "Olfactory sensory neurons" = "#FF26FF",
               
                    # muscle 
                    "Muscle progenitor cells" = "#cd5c5c",
                    "Muscle progenitor cells (Prdm1+)" = "#b22222",
                    "Myotubes" = "#8b0000",
                    # Hepatocytes
                    "Hepatocytes" = "#daa520",
                    # Lung 
                    "Lung progenitor cells" = "#ffdead",
                    # Eye and other 
                    "Pancreatic islets" = "#20B2AA",
                    "Suprachiasmatic nucleus" = "#008080",
                    # Megakaryocytes 
                    "Megakaryocytes" = "#2E8B57",
                    # Blood 
                    "Hematopoietic stem cells (Cd34+)" = "#006400",
                    "Definitive early erythroblasts (CD36-)" = "#228B22",
                    "Border-associated macrophages" = "#32CD32",
                    "Border-associated macrophages (Ms4a8a+)" = "#3CB371",

                    # CNS neurons
                    "GABAergic neurons" = "#783B3B",
                    "GABAergic cortical interneurons" = "#5C2A2A",
                    "Glutamatergic neurons" = "#DB6B6B",
                    "Neural progenitor cells (Neurod1+)" = "#E38C8C",
                    # Ependymal cells
                    "Choroid plexus" = "#6C350F",
                    # Endothelium 
                    "Brain capillary endothelial cells" = "#AF5931",
                    "Liver sinusoidal endothelial cells" = "#FFA500",
                    
                   # Intestine
                    "Extraembryonic visceral endoderm" = "#f6c100",
                    "Midgut/Hindgut epithelial cells" = "#ffecb3"
)
  
celltype_updated_colors = c(
                    ## Neural ectoderm  (brown-orange)
                    "Anterior floor plate" = "#a14a57",
                    "Floorplate and p3 domain" = "#a0522b",
                    "Anterior roof plate" = "#FFA69E",
                    "Posterior roof plate" = "#ffb6a0",
                    "Spinal cord ventral progenitors" = "#e97451",
                    "Spinal cord motor neurons" = "#ffaf87",
                    "Spinal cord dorsal progenitors" = "#e3b778",
                    "Spinal cord/r7/r8" = "#ffdcad",
                    "NMPs and spinal cord progenitors" = "#9b1c31",
                    "Dorsal telencephalon" = "#ff7f50",
                    "Telencephalon" = "#c50061",
                    "Posterior Forebrain / Diencephalon" = "red",
                    "Midbrain" = "#c1996b",
                    "Hindbrain" = "#d94e1f",
                    "Midbrain-hindbrain boundary" = "#ffae42",
                    "Hypothalamus" = "#CC3131",
                    "Hypothalamus (Sim1+)" = "#D34A4A",
                    "Anterior Forebrain" = "#ff4500",
                    "Cranial motor neurons" = "#dc143c",
                    "Neural crest (PNS neurons)" = "#732507",
                    "Neural crest (PNS glia)" = "#8b4513",

                    ## Surface ectoderm
                    "Olfactory epithelial cells" = "#ff1493",
                    "Otic epithelial cells" = "#ff69b4",
                    "Pre-epidermal keratinocytes" = "#ffb6c1",

                    # Mesoderm
                    # Paraxial mesoderm
                    "Somites / Dermomyotome" = "#97e2ff",
                    "Somites / Sclerotome" = "#156888",
                    "Mesodermal progenitors (Tbx6+)" = "#6689a1",
                    # Intermediate mesoderm
                    "Anterior intermediate mesoderm" = "#4f4a83",
                    "Lateral plate and intermediate mesoderm" = "#1d3557",
                    # Lateral plate / cardiovascular
                    "First heart field" = "#c600cb",
                    "Second heart field" = "#800080",
                    "Endocardial cells" = "#E3B5E3",
                    "Arterial endothelial cells" = "#da70d6",
                    "Endothelium" = "#d1f1d2",
                    "Pericytes" = "#90ee90",
                    "Hematoendothelial progenitors" = "#6c8e23",
                    "Primitive erythroid cells" = "#556b2f",

                    # Other mesodermal
                    "Facial mesenchyme" = "#4b0082", 
                    "FM_1" = "#6e5586",
                    "FM_0" = "#b084c0",
                    "FM_2" = "#4b0082",
                    "Chondrocytes (Atp1a2+)" = "#8a2be2",
                     "Limb mesenchyme progenitors" = "#9370db",


                    # Endoderm
                    "Gut" = "#f6f201",
                    "Ciliated nodal cells" = "#ffd700",
                    

                    # Germline
                    "Primordial germ cells" = "#ff8dfc",

                    # epithelial
                    "Notochord" = "#1242fd",
                    "Olfactory pit cells" = "#1e90ff",
                    "Granular keratinocytes" = "#87cefa",
                    "Pancreatic acinar cells" = "#00bfff",
                    "Pituitary/Pineal gland progenitors" = "#4682b4",
                    "Apical ectodermal ridge" = "#5f9ea0",
                    "Amniotic ectoderm" = "#ffa07a",
                    "Placodal area" = "#DDA0DD",
                    "Posterior intermediate mesoderm" = "#7B68EE",

                    # Neural crest PNS neurons
                    "Otic sensory neurons" = "#FF69B4",
                    "Sympathetic neurons" = "#FF47A9",
                    "Enteric neurons" = "#A30057",

                    "Olfactory sensory neurons" = "#FF26FF",
               
                    # muscle 
                    "Muscle progenitor cells" = "#cd5c5c",
                    "Muscle progenitor cells (Prdm1+)" = "#b22222",
                    "Myotubes" = "#8b0000",
                    # Hepatocytes
                    "Hepatocytes" = "#daa520",
                    # Lung 
                    "Lung progenitor cells" = "#ffdead",
                    # Eye and other 
                    "Pancreatic islets" = "#20B2AA",
                    "Suprachiasmatic nucleus" = "#008080",
                    # Megakaryocytes 
                    "Megakaryocytes" = "#2E8B57",
                    # Blood 
                    "Hematopoietic stem cells (Cd34+)" = "#006400",
                    "Definitive early erythroblasts (CD36-)" = "#228B22",
                    "Border-associated macrophages" = "#32CD32",
                    "Border-associated macrophages (Ms4a8a+)" = "#3CB371",

                    # CNS neurons
                    "GABAergic neurons" = "#783B3B",
                    "GABAergic cortical interneurons" = "#5C2A2A",
                    "Glutamatergic neurons" = "#DB6B6B",
                    "Neural progenitor cells (Neurod1+)" = "#E38C8C",
                    # Ependymal cells
                    "Choroid plexus" = "#6C350F",
                    # Endothelium 
                    "Brain capillary endothelial cells" = "#AF5931",
                    "Liver sinusoidal endothelial cells" = "#FFA500",
                    
                   # Intestine
                    "Extraembryonic visceral endoderm" = "#f6c100",
                    "Midgut/Hindgut epithelial cells" = "#ffecb3"
)


day_colors = c("E8.0-E8.5" = "#78c1e2",
               "E8.75" = "#2c7bb6",
               "E9.0" = "#75b38d",
               "E9.25" = "#2e5d37",
               "E9.5" = "#c786cd",
               "E9.75" =  "#8f3b96"
)

celltype_order <- c(
  "Telencephalon" = "Telencephalon (1)",
  "Posterior Forebrain / Diencephalon" = "Posterior Forebrain / Diencephalon (2)",
  "Midbrain" = "Midbrain (3)",
  "Midbrain-hindbrain boundary" = "Midbrain-hindbrain boundary (4)",
  "Hindbrain" = "Hindbrain (5)",
  "Floorplate and p3 domain" = "Floorplate and p3 domain (6)",
  "Spinal cord/r7/r8" = "Spinal cord/r7/r8 (7)",
  "NMPs and spinal cord progenitors" = "NMPs and spinal cord progenitors (8)",
  "Notochord" = "Notochord (9)",
  "Mesodermal progenitors (Tbx6+)" = "Mesodermal progenitors (Tbx6+) (10)",
  "Somites / Sclerotome" = "Sclerotome (11)",
  "Somites / Dermomyotome" = "Dermomyotome (12)",
  "Gut" = "Gut (13)",
  "Placodal area" = "Placodal area (14)"
)

cell_type_colored_numbered <- setNames(
  cell_type_colored[names(celltype_order)],  # get the colors for those names (in correct order)
  celltype_order                             # new names, e.g. "Telencephalon (2)"
)
