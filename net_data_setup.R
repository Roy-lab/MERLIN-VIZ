### App Setup! 

library(tidyverse)
library(tidygraph)
library(pracma)
library(DT)
library(Matrix)
library(data.table)

source('aux_functions.R')

### Files used for netData generation. 
prefix <- "/Volumes/"


all_nodes_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2026/suvo_work/spencer_style_network/unique_nodes_v2.txt")
edge_list_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2025/marina_work/results/Merlinp_results/Lambda_0100/consensus/n20_subsamples_lambda_0100_0_8.txt")
module2gene_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2025/marina_work/results/Merlinp_results/Lambda_0100/consensus/consensus_module_0_2_geneset.txt")
module_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2025/marina_work/results/Merlinp_results/Lambda_0100/consensus/consensus_module_0_2_geneset_enrichAnalyzer.txt")
regulator_enrich_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2025/marina_work/results/Merlinp_results/Lambda_0100/consensus/regulator_enrichAnalysis_0_2_details.txt")
go_enrich_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2025/marina_work/results/Merlinp_results/Lambda_0100/consensus/go_enrichAnalysis_0_2_details.txt")
gene2genename_file <- NULL
Ortholog_1_to_1_file <- NULL
Ortholog_file <- NULL
go_file = paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2025/marina_work/data/mousegotermap_regnet.txt")
gene_desc_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2026/suvo_work/spencer_style_network/consensus_module_0_2_geneset_names.txt")
regulator_list_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2026/suvo_work/spencer_style_network/net1_transcription_factors.txt")
expression_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2025/marina_work/results/Nca/Lambda_0100/Merlinp_inputs/net1_expression_with_header_gene_by_cell.txt") #header "Gene"
cell_mapping_file <- paste0(prefix, "/wid/projects7/Roy-singlecell2/bookchapter_MERLIN_2026/suvo_work/spencer_style_network/sample_annotation.txt")

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
                          grouping_indices = grouping_indice, 
                          outfile = 'test.Rdata')


