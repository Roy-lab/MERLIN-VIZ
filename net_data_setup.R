### App Setup! 

library(tidyverse)
library(tidygraph)
library(pracma)
library(DT)
library(Matrix)
library(data.table)

source('aux_function.R')

### Files used for netData generation. 
prefix <- "/Volumes/"
all_nodes_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/js_work/Saha_u01/Paper_analysis_figures_2025/merlin_shiny_app/PR_UTF_2D3D/all_genes.txt")
edge_list_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/sahalab/results/UTF_2D3D_merlin/PR_15_17_20/output/processed/nets/net_0.6.txt")
module2gene_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/sahalab/results/UTF_2D3D_merlin/PR_15_17_20/output/processed/modules/pseudobulk_modules/module.0.3_geneset_zeromean_reorder.txt")
module_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/js_work/Saha_u01/Paper_analysis_figures_2025/merlin_visualization/UTF_PR_merlin/gnset_maxdim.txt")
go_file = paste0(prefix, "/wid/projects2/Roy-common/data/data_new/human/go/hg38/hg38_goterms_regnet.txt")
regulator_enrich_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/sahalab/results/UTF_2D3D_merlin/PR_15_17_20/output/processed/modules/enrichments/merlin.0.6.0.3_details.txt")
go_enrich_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/js_work/Saha_u01/Paper_analysis_figures_2025/merlin_visualization/UTF_PR_merlin/gnset_gobp_enrich_details.txt")
gene2genename_file <- NULL
Ortholog_1_to_1_file <- NULL
Ortholog_file <- NULL
gene_desc_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/data/gene_annotations/gene_names.txt")
regulator_list_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/sahalab/ref/human_regulators_nca.txt")
expression_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/sahalab/results/UTF_2D3D_merlin/PR_15_17_20/output/processed/modules/pseudobulk_modules/UTF_2D3D_PR_15_17_20_CxGsubmat_zeromean.txt")
cell_mapping_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/js_work/Saha_u01/Paper_analysis_figures_2025/merlin_visualization/utf_sample_reference.txt")



################### Make R Data Files *Only need to run once###################
expression_data <- prepareExpression(expression_file, save_struct = FALSE)
grouping_indice <- prepareCellMapping(expression_data, cell_mapping_file, save_struct = FALSE)

makePostProcessDataStruct(all_nodes_file, 
                          edge_list_file,
                          module2gene_file = module2gene_file, 
                          go_file = go_file,
                          module_file = module_file, 
                          regulator_enrich_file = regulator_enrich_file,
                          go_enrich_file = go_enrich_file, 
                          gene_desc_file = gene_desc_file, 
                          regulator_list_file = regulator_list_file, 
                          expression_data = expression_data, 
                          grouping_indices = grouping_indice)


