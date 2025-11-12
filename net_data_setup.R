### App Setup! 

library(tidyverse)
library(tidygraph)
library(pracma)
library(DT)
library(Matrix)
library(data.table)

source('aux_functions.R')

### Files used for netData generation. 
prefix <- "/mnt/dv/"
all_nodes_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/expression_by_cluster/union_genes_no_cl_0_3.txt")
edge_list_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/Merlin_results/consensus/network_0_8.txt")
module2gene_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/Merlin_results/consensus/consensus_module_0_3_geneset.txt")
module_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/Merlin_results/consensus/consensus_module_0_3_geneset_enrichAnalyzer.txt")
go_file = paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/data/Changa_sorghum/sorghum_scRNAseq_data_to_share_Roy_team/GO_regnet_v2.txt")
regulator_enrich_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/Merlin_results/consensus/regulator_enrichAnalysis_0_3_details.txt")
go_enrich_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/Merlin_results/consensus/GO_enrichAnalysis_0_3_details.txt")

gene2genename_file <- NULL
Ortholog_1_to_1_file <- NULL
Ortholog_file <- NULL

gene_desc_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/gene_names.txt")
regulator_list_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/Merlin_inputs/union_regulators_intersected.txt")
expression_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/expression_by_cluster/filtered_merged_expr.txt") #header "Gene"
cell_mapping_file <- paste0(prefix, "/wid/projects7/Roy-plants/kirstlab/results/Changa_MERLIN/cluster_info_for_expr.txt")



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


