### App Setup! 

library(tidyverse)
library(tidygraph)
library(pracma)
library(DT)
library(Matrix)
library(data.table)

source('aux_functions.R')

### Files used for netData generation.
# Use example_data_minimal (git-friendly subset: Egr1 hub + Cluster1002 module)
example_dir <- "example_data_minimal"

all_nodes_file        <- file.path(example_dir, "unique_nodes_v2.txt")
edge_list_file        <- file.path(example_dir, "n20_subsamples_lambda_0100_0_8.txt")
module2gene_file      <- file.path(example_dir, "consensus_module_0_2_geneset.txt")
module_file           <- file.path(example_dir, "consensus_module_0_2_geneset_enrichAnalyzer.txt")
regulator_enrich_file <- file.path(example_dir, "regulator_enrichAnalysis_0_2_details.txt")
go_enrich_file        <- file.path(example_dir, "go_enrichAnalysis_0_2_details.txt")
gene2genename_file    <- NULL
Ortholog_1_to_1_file  <- NULL
Ortholog_file         <- NULL
go_file               <- file.path(example_dir, "mousegotermap_regnet.txt")
gene_desc_file        <- file.path(example_dir, "consensus_module_0_2_geneset_names.txt")
regulator_list_file   <- file.path(example_dir, "net1_transcription_factors.txt")
expression_file       <- file.path(example_dir, "net1_expression_with_header_gene_by_cell.txt") #header "Gene"
cell_mapping_file     <- file.path(example_dir, "sample_annotation.txt")

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
                          outfile = 'net_data.Rdata')


