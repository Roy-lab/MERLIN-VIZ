library(tidyverse)
library(tidygraph)
library(pracma)
library(DT)


### Files used for netData generation. 
#prefix <- "/mnt/dv/"
#all_nodes_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/data/katie_t_cells/genes_filtered.txt", sep = "")
#edge_list_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/results/MERLIN/katie_t_cells/consensus/network_0.8.txt", sep = "")
#module2gene_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/results/MERLIN/katie_t_cells/consensus/consensus_module_0.1_geneset.txt", sep = "")
#module_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/results/MERLIN/katie_t_cells/consensus/consensus_module_0.1_geneset_enrichAnalyzer.txt", sep = "")
#go_file = paste(prefix, "wid/projects2/Roy-common/data/data_new/human/go/hg38/hg38_goterms_regnet.txt", sep = "")
#regulator_enrich_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/results/MERLIN/katie_t_cells/consensus/regulator_enrichAnalysis_0.1_details.txt", sep = "")
#go_enrich_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/results/MERLIN/katie_t_cells/consensus/enrichAnalysis_0.1_details.txt", sep = "")
#gene2genename_file <- NA
#gene_desc_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/data/gene_annotations/gene_names.txt", sep = "")
#regulator_list_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/data/katie_t_cells/regulators_intersected_with_transgene.txt", sep = "")
#expression_file <- paste(prefix, "wid/projects7/Roy-singlecell2/sahalab/Tcell/jeremy_work/data/katie_t_cells/expression_geneids_subset_filtered_barcodes_transposed.txt", sep = "")

################### MERLIN_VIZ DEFAULTS FOR SCALES #############################
title <- "Katie T-cell"
default_edge_color_pallette <- "RdBu"
default_node_color_pallette <- "Pastel1"
default_gene <- "AACS"
default_expression_heatmap <- "Reds"
default_expression_range <- c(0, 5)
default_expression_min <- 0
default_expression_max <- 10
default_tfa_palette_heatamp <- "PiYG"
default_tfa_range <- c(-2, 2)
default_tfa_min <- -10
default_tfa_max <- 10

################### Make R Data Files *Only need to run once###################
makePostProcessDataStruct <- function (all_nodes_file, edge_list_file,
                          module2gene_file, go_file, module_file, 
                          regulator_enrich_file, go_enrich_file, 
                          gene2genename_file, gene_desc_file, regulator_list_file, expression_file, max_samps = 5000)
{
###### Generate labeled Node Set
genes2modules <- read_tsv(module2gene_file, c("feature", "module"))
genes2modules <- genes2modules %>% 
  group_by(module) %>% 
  summarize(count = n()) %>% 
  right_join(genes2modules, by = "module") %>%
  ungroup() %>% 
  filter(count > 4) %>%
  select(feature, module)

go <- read_tsv(go_file, col_names = c("feature", "go"), skip = 1) 

go_IC <- go %>% 
  group_by(go) %>%
  summarise(freq = n()) %>%
  mutate(IC = -log(freq/sum(freq)))

go <- left_join(go, go_IC) %>%
  group_by(feature) %>%
  arrange(desc(IC), .by_group = TRUE) %>% 
  select(feature, go) %>% 
  chop(go) %>% 
  ungroup() %>% 
  arrange(feature)

regulators <- read_tsv(regulator_list_file, col_names = FALSE) %>%
  rename("feature" = "X1") %>%
  mutate(regulator = TRUE)



if(!is.na(gene2genename_file)){
  gene_map <- read_tsv(gene2genename_file, col_names = FALSE) %>%
    rename("feature" = "X1", "Common Name" = "X2")
}

gene_desc <- read_tsv(gene_desc_file, col_names = c("feature", "Description")) #%>%
#rename("feature" = "X1", "Description" = "X2")

nodes <- read_tsv(file = all_nodes_file, col_names = "feature") %>% 
  rowid_to_column("id") %>%
  left_join(genes2modules) %>%
  left_join(go) %>%
  #left_join(gene_map) %>%
  left_join(gene_desc) %>%
  left_join(regulators) 

if(is.na(gene2genename_file))
{
  nodes <- nodes %>% 
    mutate(`Common Name` = feature) 
}else
{
  nodes <- nodes %>% 
    left_join(gene_map, by = "feature") %>%
    mutate(`Common Name` = ifelse(is.na(`Common Name`), feature, `Common Name`))
  
}

nodes$regulator <- sapply(nodes$regulator, function(x){ifelse(is.na(x), 'tar', 'scr')})
nodes$module <- sapply(nodes$module, function(x){ifelse(is.na(x), -9999, x)})


nca_idx <- which(grepl('_nca', nodes$feature))
for(idx in nca_idx){
  node_name <- nodes$feature[idx]
  str_name <- str_remove(node_name, '_nca')
  str_idx <- which(str_name == nodes$feature)
  nodes$go[idx] <- nodes$go[str_idx]
  nodes$`Common Name`[idx] <- str_c(nodes$`Common Name`[str_idx], '_nca')
  nodes$Description[idx] <-nodes$Description[str_idx]
} 

gene_map <- nodes %>% filter( !grepl('AFUA_', `Common Name`)) 
common_name = gene_map$`Common Name`
feature_name = gene_map$feature
genename_map <- tibble(common_name, feature_name)

## load in Expression matrix 
expression <- vroom::vroom(expression_file) %>% 
  select(Gene, where(~ !all(is.na(.x) | .x == ""))) %>%
  pivot_longer(cols = !Gene) %>% 
  group_by(Gene) %>%
  dplyr::summarize(expression = list(value))


### Generate Edge Set
routes <- read_tsv(edge_list_file, c("source", "target", "weight"))
edges <- routes %>% 
  left_join(nodes, by = c("source" = "feature")) %>% 
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("target" = "feature")) %>% 
  rename(to = id)

edges <- edges %>%
  select('from', 'to', 'weight')

### Make Tidy Graph Structure 
Net <- tbl_graph(nodes= nodes, edges = edges)
Net <- Net %>%
  #convert(to_undirected) %>%
  mutate(neighbors = map_local(order = 1, .f = function(neighborhood, node, ...) {
         as_tibble(neighborhood, active = 'nodes')$feature
      }))
degree_v <-  Net %N>% as_tibble() %>% rowwise() %>% summarize(length(neighbors)) -1
Net <- Net %>% mutate(degree = degree_v$`length(neighbors)`)

Net <- left_join(Net, expression, by = c("feature" = "Gene"))

## add correlation between expression vectors for all edges
Net <- Net %E>% 
  mutate(Correlation = map2_dbl( 
    .N()$expression[from],
    .N()$expression[to], 
     ~ cor(.x, .y))
    ) 

safe_slope <- function(x, y) {
  v <- var(x)
  if (v == 0 || is.na(v)) return(NA_real_)
  cov(x, y) / v
}
## Add regression weights betwee all edges.   
Net <- Net %E>%
  mutate(Reg_weight = map2_dbl(
    .N()$expression[from], 
    .N()$expression[to], 
    safe_slope
  )) 


## Subsample expression after for visualization only: 

set.seed(42)  
sample_size <- max_samps
full_length <- length((Net %N>% pull(expression))[[1]])
subsample_indices <- round(seq(1, full_length, length.out = sample_size))

Net <- Net %>%
  activate(nodes) %>%
  mutate(expression = map(expression, ~ .x[subsample_indices]))

### Generate Module Structure 
Module <- read_tsv(module_file, c("module", "gene_list"))  
Module$gene_list <- str_split(Module$gene_list, '#')
go_enrich <- read_tsv(go_enrich_file, col_names = c("module", "go", 
                      "go_pvalue", "go_correct_pvalue", "go_num_tot_genes", 
                      "num_tot_go_genes", "num_module_genes", 
                      "num_module_go_genes", "go_enrichment", "go_genes"))
go_enrich$go_genes <- str_split(go_enrich$go_genes, ';') 
go_enrich <- left_join(go_enrich, go_IC) %>%
  group_by(module) %>% arrange(desc(IC), .by_group = TRUE) %>% 
  mutate(IC = NULL, freq = NULL) %>% 
  nest(GO=c("go", "go_pvalue", "go_correct_pvalue", 
                                  "go_num_tot_genes", "num_tot_go_genes", 
                                  "num_module_genes", "num_module_go_genes", 
                                  "go_enrichment", "go_genes"))

merlin_enrich <- read_tsv(regulator_enrich_file, col_names = c("module", 
                          "regulator", "reg_pvalue", "reg_correct_pvalue", 
                          "reg_num_tot_genes", "num_tot_regulator_genes", 
                          "reg_num_module_genes", "num_module_regulator_genes", 
                          "regulator_enrichment", "reg_target_genes"))
regulator2module_map <- merlin_enrich %>% select("module", "regulator") %>% mutate(module = strtoi(str_remove(module, 'Cluster'))) %>% chop(module) %>% rename("enriched_modules" = module)
Net <- Net %N>% left_join(regulator2module_map, by = c("feature" = "regulator"))
merlin_enrich$reg_target_genes <- str_split(merlin_enrich$reg_target_genes, ';')
merlin_enrich <- nest(merlin_enrich, "regulators"= c("regulator", "reg_pvalue", 
                          "reg_correct_pvalue", "reg_num_tot_genes", 
                          "num_tot_regulator_genes", "reg_num_module_genes", 
                          "num_module_regulator_genes", "regulator_enrichment", 
                          "reg_target_genes") )

Module <- Module %>% 
  left_join(go_enrich, by = "module") %>%
  left_join(merlin_enrich, by = "module")

Module$module <- strtoi(str_replace_all(Module$module, "Cluster", ""))

module_ids <- Module %>% pull(module)
enrich_2_module <- Module %>% select(module, GO) %>% unnest(GO) %>% select(go, module) %>% group_by(go) %>% chop(module)

enriched_go_terms <- enrich_2_module %>%
pull(go)

 
genes <- unique(Net %N>% as_tibble() %>% pull(feature))

save(list = c("Net", "Module", "enriched_go_terms", "module_ids","enrich_2_module", "genes", "genename_map"), file = "net_data.Rdata")

return(list(Net, Module, enriched_go_terms, module_ids, enrich_2_module, genes, genename_map)) 
}


makeLaplacian <- function(Net){
  degree_vec <- Net %N>% pull(degree)
  D <- diag(degree_vec)
  adj = matrix(0, length(degree_vec), length(degree_vec))
  Edges <- Net %E>% as_tibble()
  from <- Edges %>% pull(from)
  to <- Edges %>% pull(to)
  for(i in 1:length(from)){
    adj[from[i], to[i]]= 1
    adj[to[i], from[i]]= 1
  }
  L <- D - adj;
  save(list = c("Net", "Module", "enriched_go_terms", "module_ids","enrich_2_module", "genes", "L"), file = "net_data.Rdata")
  return(L)
}


MakeKernel <- function(L, lambda){
  num_nodes <- size(L)[1]
  I <- eye(num_nodes)
  inside <- I + lambda*L
  kernel <- inv(inside)
  return(kernel)
}

###############################################################################
##################### Search Functions ########################################
###############################################################################
searchForModule<- function(Module, moduleID){
  modIdx<-which(Module$module == moduleID)
  gene_in_module <- Module$gene_list[[modIdx]]
  regulators <- Module$regulators[[modIdx]]$regulator
  mod_genes<- c(gene_in_module, regulators)
return(mod_genes)
}

searchForGene <-function(Net,Module, gene){
  moduleIDs <- Net %N>%
    filter(feature %in% gene) %>%
    pull(module)
  moduleIDs <- unique(moduleIDs[!is.na(moduleIDs)])
  genes <-list()
  for (ID in moduleIDs) {
    new_genes <- searchForModule(Module, ID )
    genes <- append(unlist(genes), unlist(new_genes))
  }
  
  neighbors <- Net %N>%
    filter(feature %in% gene) %>%
    pull(neighbors)
  
  genes <- append(unlist(genes), unlist(unlist(neighbors)))
  print(unique(genes))
  return(unique(genes))
}

searchForGeneList <-function(Net, Module, gene_list, search_additional){
  result_list <- Net %N>%
    filter(feature %in% gene_list) %>%
    pull(feature)
  
  if("mod" %in% search_additional){
    mod_genes <-list()
    moduleIDs <- unique(c(
      Net %N>%
        filter( feature %in% result_list) %>%
        filter( module != -9999) %>%
        pull(module), 
      unlist(Net %N>%
        filter( feature %in% result_list) %>%
        pull(enriched_modules))))
    moduleIDs <- moduleIDs[!is.na(moduleIDs)]
    
    for (ID in moduleIDs) {
      new_genes <- searchForModule(Module, ID )
      mod_genes <- append(unlist(mod_genes), unlist(new_genes))
    }
    result_list <- append(unlist(result_list), unlist(mod_genes))
  }
  
  if("neigh" %in% search_additional){
    neighbors <- Net %N>%
      filter(feature %in% result_list) %>%
      pull(neighbors)
    result_list <- append(unlist(result_list), unlist(unlist(neighbors)))
  }
return(unique(result_list))
}

computeEnrichment <- function(Module, gl, num_genes){
    Module <- Module %>% rowwise() %>% 
      mutate(m_size = length(gene_list)) %>%
      mutate(intersect_size = length(intersect(gl, gene_list))) %>%
      mutate("Genes on List" = ifelse(intersect_size > 0, list(intersect(gl, gene_list)), NA)) %>% #ifelse(intersect_size > 0, intersect(gl, gene_list), NA)) %>%
      mutate(enrich_pval = phyper(intersect_size, length(gl), num_genes, m_size, lower.tail=FALSE)) %>%
      mutate(corrected_pval = ifelse(enrich_pval * dim(Module)[1] < 1, enrich_pval * dim(Module)[1], 1))
    return(Module)
}

########################## Node Diffusion #####################################

generateScoreVector <- function(Net, gene_list){
  nodes <- Net %N>% as_tibble() %>% select(feature)
  scores <- left_join(nodes, Net %N>% as_tibble() %>% 
    filter(str_detect(feature, paste(unlist(gene_list), collapse="|"))) %>%
    select(feature) %>% 
    mutate(score = 100)) %>% 
    mutate (score = replace_na(score, 0))
  return(scores$score)
}


loadScoreVector <- function(Net, score_data) {
  nodes <- Net %N>% as_tibble() %>% mutate(score = 0)
  for( i in 1:length(score_data[[1]])){ 
    nodes <- nodes %>%
      mutate(score = replace(score, 
             which(str_detect(feature, score_data %>% slice(i) %>% pull(feature))),
             score_data %>% slice(i) %>% pull(score)))
  }
  return(nodes$score)
}

computeDiffusionScore <- function (Net, score_data, kernel){
  score <- loadScoreVector(Net, score_data)
  diff_score <- kernel %*% score
  Net <- Net %N>% mutate("score" = as.vector(diff_score))
  return(Net)
}



###############################################################################
##################### Subgraph Functions ######################################
###############################################################################
moduleSubgraph <- function(Net, Module, module_id){
  mod_genes <- searchForModule(Module, module_id)
  sub_graph <-induceSubraph(Net, mod_genes) %E>%
    mutate(color_code = "#666")
  return(sub_graph)
}

geneSubgraph <- function(Net, Module, gene){
  list_genes <- searchForGene(Net, Module, gene)
  sub_graph <- induceSubraph(Net, list_genes) %E>%
    mutate(color_code = "#666")
  return(sub_graph)
}

geneListSubgraph <- function(Net, Module, gene_list, search_additional){
  if("stein" %in% search_additional & length(gene_list) > 1){
    st<-buildSteinerTrees(Net, gene_list) %E>%
      mutate(is_steiner = TRUE)
      gene_list <- st %N>% pull(feature)
    
  }
  list_genes <- searchForGeneList(Net,Module, gene_list, search_additional)
  sub_graph <- induceSubraph(Net, list_genes)
  
  if("stein" %in% search_additional & length(gene_list) > 1){
    sub_graph <- graph_join(sub_graph, st) %E>%
      mutate(is_steiner = replace_na(is_steiner, FALSE)) %>%
      mutate(color_code = if_else(is_steiner, "#fb8072", "#666"))
  }else{
    sub_graph <- sub_graph %E>% 
      mutate(color_code = "#666")
    }
  return(sub_graph)
}

goSubgraph <- function(Net, Module, enrich_2_module, go_term){
  modules_list<-unlist(enrich_2_module %>%
                         filter(go == go_term) %>%
                         pull(module))
  mod_genes = list()
  for(i in 1:length(modules_list)){
    mod_genes<- append(unlist(mod_genes), unlist(searchForModule(Module, modules_list[i])))
  }
  sub_graph <- induceSubraph(Net, mod_genes) %E>%
    mutate(color_code = "#666")
  return(sub_graph)
}

diffScoreSubgraph <- function(Net, min_targets, top_regs){
  nodes <- Net %N>% as_tibble()
  nodes <- nodes %>% arrange(desc(score)) %>% filter(regulator == 'scr') %>% filter(degree >= min_targets)
  top5 <- nodes  %>% slice(1:top_regs)
  genes <- unlist(top5$neighbors) 
  sub_graph <- induceSubraph(Net, genes) %E>%
    mutate(color_code = "#666")
  return(sub_graph)
}

induceSubraph <- function(Net, list){
  sub_graph <- Net %N>%
    convert(to_subgraph, feature %in% list)
  return(sub_graph)
}

graph2NodeEdgeTables <- function(Net){
  graph_nodes <- Net %N>%
    as_tibble() %>%
    mutate(id = id -1)
  graph_edges <- Net %E>%
    as_tibble() %>%
    mutate(from =from - 1, to = to - 1)
  return( list(graph_nodes, graph_edges))
}


###############################################################################
##################### Steiner Tree Construction ###############################
###############################################################################
getDistMatrix <- function (Net, gene_list){
gene_id <- Net %>% 
  filter(feature %in% gene_list) %>%
  pull(id)

gene_name <-  Net %>% 
	filter(feature %in% gene_list) %>%
	pull(feature)

dist_matrix <- Net %N>%
  as_tibble() %>%
  select(feature)

for(idx in 1:length(gene_id)){
  root <- gene_id[idx]
  name <- gene_name[idx]
  dist_matrix <- Net %N>% 
  mutate(dist = bfs_dist(root)) %>%
  as_tibble() %>%
  select(feature, dist) %>%
  rename(!!name := dist) %>%
  right_join(dist_matrix, by="feature")  
}

dist_matrix <- dist_matrix %>% replace(. ==0 , length(dist_matrix$feature)+1)
return(dist_matrix)
}


buildSteinerTrees <- function(Net, gene_list){
  Net <- Net %>% convert(to_undirected)
  dist_matrix <- getDistMatrix(Net, gene_list)
  gene_names <- colnames(dist_matrix)[2:length(dist_matrix)]
  dist2graph <- tibble(gene_names)
  dist2graph <- dist2graph %>%
    mutate(Dist = Inf) %>%
    mutate(Closest = "")
  
  
## Find Closest Nodes
  steiner_tree=tbl_graph()
  restrict_dist_matrix <- dist_matrix %>%
    filter(feature %in% gene_names) 
  
  for(i in 2:length(restrict_dist_matrix)){
    gene <- colnames(restrict_dist_matrix[i])
    idx <-which(dist2graph$gene_names == gene)
    if(all(is.na(restrict_dist_matrix[i]))){
      dist2graph$Closest[idx[1]] <- gene
      dist2graph$Dist[idx[1]] <- length(dist_matrix$feature)+1
    }else{
      dist <- min(restrict_dist_matrix[i], na.rm = TRUE)
      match_idx <- which(restrict_dist_matrix[i] == min(restrict_dist_matrix[i], na.rm = TRUE))
      closest <- restrict_dist_matrix$feature[match_idx[1]]
      dist2graph$Dist[idx[1]] <- dist
      dist2graph$Closest[idx[1]] <- closest
    }
  }
  
## Select minimum length path   
  path_2_add<-dist2graph %>% 
    arrange(Dist) %>%
    slice(1) %>%
    select("gene_names", "Closest")
  path_2_add<-c(path_2_add$gene_names,path_2_add$Closest)
  
## Find id of node_2_add  
  nodes_2_add<- Net %N>% 
    filter(feature %in% path_2_add) %>%
    pull(id)
    

## Find path and initialize Steiner tree. 
  steiner_tree <- Net %N>%
    convert(to_shortest_path, nodes_2_add[1], nodes_2_add[2])
  
  nodes <- steiner_tree %N>%
    pull(feature)
  
## Prune search matrix & dist2graph    
  dist_matrix <- dist_matrix %>% 
    select(-intersect(gene_list, nodes))
  dist2graph <- dist2graph %>%
    filter(gene_names %in% setdiff(gene_names, nodes))
  
    
## Now we to find distance to nodes in graph.
  while(length(intersect(nodes, gene_names)) < length(gene_names)){
    restrict_dist_matrix <- dist_matrix %>%
      filter(feature %in% nodes)

    for(i in 2:length(restrict_dist_matrix)){
      gene <- colnames(restrict_dist_matrix[i])
      idx <-which(dist2graph$gene_names == gene)
      if(all(is.na(restrict_dist_matrix[i]))){
        dist2graph$Closest[idx[1]] <- gene
        dist2graph$Dist[idx[1]] <- length(dist_matrix$feature)+1
      }else{
        dist <- min(restrict_dist_matrix[i], na.rm = TRUE)
        match_idx <- which(restrict_dist_matrix[i] == min(restrict_dist_matrix[i], na.rm = TRUE))
        closest <- restrict_dist_matrix$feature[match_idx[1]]
        dist2graph$Dist[idx[1]] <- dist
        dist2graph$Closest[idx[1]] <- closest
      }
    }
  
    path_2_add<-dist2graph %>% 
      arrange(Dist) %>%
      slice(1) %>%
      select("gene_names", "Closest")
    path_2_add<-c(path_2_add$gene_names,path_2_add$Closest)
  
    nodes_2_add<- Net %N>% 
      filter(feature %in% path_2_add) %>%
      pull(id)
    
    if( length(nodes_2_add) == 1){
      node_info<-Net %N>%
        filter(id %in% nodes_2_add) %>%
        as_tibble()
      steiner_tree<- steiner_tree %>%
        bind_nodes(node_info)
    }else{
      Path <- Net %N>%
        convert(to_shortest_path, nodes_2_add[1], nodes_2_add[2])
       
        
    
      steiner_tree <- steiner_tree %>% 
        graph_join(Path)
    }
    
    nodes <- steiner_tree %N>%
      pull(feature)
    
    ## Prune search matrix & dist2graph    
    dist_matrix <- dist_matrix %>% 
      select(setdiff(colnames(dist_matrix), nodes))
    dist2graph <- dist2graph %>%
      filter(gene_names %in% setdiff(gene_names, nodes))
  }
return(steiner_tree)
}

###############################################################################
######################### Print Functions #####################################
###############################################################################


printNodeInfo <- function(Net, node_name){
  if(is.na(node_name)){
    return(" ")
  }else{
  node_data <- Net %N>%
    filter(feature == node_name) %>%
    as_tibble()
  node_id <- node_data %>% pull(id)
  
  neighbors <- setdiff(Net %>%
                          convert(to_local_neighborhood, node_id) %N>%
                          as_tibble() %>%
                          pull(feature), node_name)
  
  text <- sprintf("Node Name: %s<br/>Module: %d<br/>GO Terms:  %s<br/>Neighbors: %s", 
                  node_name, node_data$module, 
                  paste(unlist(node_data$go), collapse = ', '), 
                  paste(unlist(neighbors), collapse = ', '))  
  return(text)
  }
}

printModuleInfo <- function(Module, module_id, gene_list, genes){
  if(is.na(module_id)){
    return(" ")
  }else{
    
    if(!is_empty(gene_list)){
      Module<-computeEnrichment(Module, gene_list, length(genes))
    }
    
    module_info <- Module %>%
      filter(module == module_id)
    
    module_regulators <- module_info %>%
      select(regulators) %>%
      unnest(regulators)
    if(nrow(module_regulators) == 0){
      regulator_text <- "No enriched regulators"
    }else{
      regulator<- module_regulators$regulator[1]
      p_value <- module_regulators$reg_correct_pvalue[1]
      target <- module_regulators$reg_target_genes[[1]]
      regulator_text <-sprintf ("Regulator Name: %s &emsp;P-Value: %.02e<br/> Target Genes: %s<br/>", regulator, p_value, paste(unlist(target), collapse = ', '))
      if(length(module_regulators$regulator) > 1){
        for(i in 2:length(module_regulators$regulator)){
          regulator<- module_regulators$regulator[i]
          p_value <- module_regulators$reg_correct_pvalue[i]
          target <- module_regulators$reg_target_genes[[i]]
          regulator_text <-sprintf ("%sRegulator Name: %s &emsp;P-Value: %.02e<br/> Target Genes: %s<br/>", regulator_text, regulator, p_value, paste(unlist(target), collapse = ', '))
        }
      }
    }
    
    module_go <- module_info %>% 
      select(GO) %>%
      unnest(GO) 
    
    if(nrow(module_go) == 0){
      go_text <- "No enriched GO terms"
    }else{
      go <- module_go$go[1]
      p_value <- module_go$go_correct_pvalue[1]
      go_genes <- module_go$go_genes[[1]]
      go_text <-sprintf ("Go Term: %s &emsp;P-Value: %.02e<br/> GO Genes: %s<br/>", go, p_value, paste(unlist(go_genes), collapse = ', '))
      if(length(module_go$go) >1 ){
        for(i in 2:length(module_go$go)){
          go <- module_go$go[i]
          p_value <- module_go$go_correct_pvalue[i]
          go_genes <- module_go$go_genes[[i]]
          go_text <-sprintf ("%sGo Term: %s &emsp;P-Value: %.02e<br/> GO Genes: %s<br/>", go_text, go, p_value, paste(unlist(go_genes), collapse = ', '))
        }
      }
    }
    
    if(is_empty(gene_list)){
      text <- sprintf("Module Name: %s<br/>Module Genes: %s<br/><br/>Enriched Regulators:<br/>%s<br/><br/>Enriched GO:<br/>%s <br/><br/>", 
                      module_id, paste(unlist(module_info$gene_list), collapse = ', '), 
                      regulator_text, 
                      go_text)  
      
    }else{
      text <- sprintf("Module Name: %s<br/>module enrichment p-value: %.02e<br/>module corrected p-value: %.02e<br/>Genes from gene list: %s<br/>Module Genes: %s<br/><br/>Enriched Regulators:<br/>%s<br/><br/>Enriched GO: %s <br/><br/>", 
                      module_id, 
                      module_info$enrich_pval, 
                      module_info$corrected_pval, 
                      paste(ifelse( length(intersect(module_info$gene_list[[1]], gene_list)) > 0, unlist(intersect(module_info$gene_list[[1]], gene_list)), ""), collapse = ', '),
                      paste(unlist(module_info$gene_list), collapse = ', '),
                      regulator_text, 
                      go_text)  
    }
    return(text)
  }
}

printAllModuleInfo <- function(SubNet, Module, gene_list, genes){
 if(!is_empty(gene_list)){
   text_info<-""
   unique_modules <-unique(SubNet %N>% 
                             as_tibble() %>%
                             pull(module))
   
   for(id in unique_modules){
     text_info<- sprintf('%s %s', text_info, printModuleInfo(Module, id, gene_list, genes))
   }
 }else{
  text_info<-""
  unique_modules <-unique(SubNet %N>% 
   as_tibble() %>%
   pull(module))
  for(id in unique_modules){
   text_info<- sprintf('%s %s', text_info, printModuleInfo(Module, id, list()))
  }
 }
 return(text_info)
}

getModuleID <- function(Net, node_name){
  module_id <- Net %N>%
    filter(feature == node_name) %>%
    as_tibble() %>% 
    pull(module)
  print(module_id)
  return(module_id)
}

prepNodeTable <- function(Nodes_Table, disp_num){
  Nodes_Table <- Nodes_Table %>% mutate(.tidygraph_node_index = NULL, enriched_modules = NULL) %>%
    rowwise() %>% 
    mutate(go = paste(unlist(setdiff(go[1:disp_num], NA)), collapse = ' <br/>' )) %>% 
    mutate(neighbors = paste(unlist(sapply(neighbors, function(x) 
    if(x %in% genename_map$feature_name){
      x <- str_pad(genename_map$common_name[which(x == genename_map$feature_name)], 12, side = "both")
    }else{
      x <- x
    })), collapse = ' | ')) %>% 
    #select(!geneSuper) %>% 
    select(!expression) %>%
    rename("Gene Name" = "feature") %>%
    mutate("Gene Name" = sprintf('<a href=https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s target="_blank" rel="noopener noreferrer"> %s</a>', str_replace(`Gene Name`, '_nca', ''), `Gene Name`)) %>%
    mutate("id" = NULL) %>%
    mutate("regulator" = NULL)
  
  Nodes_Table$module[which(Nodes_Table$module == -9999)] <- NA
  return(Nodes_Table)
}

prepModuleTable <- function(Module_Table, method, disp_num  = 5 ){
  GO <- Module_Table %>% select(module, GO) %>%
    unnest(GO)
  GO <- GO %>% group_by(module) %>% slice_head(n = disp_num)
  
  regulators <- Module_Table %>% select(module, regulators) %>%
    unnest(regulators)
  
  GO <- GO %>% rowwise%>% 
    mutate(GO = ifelse(length(module) > 0, sprintf("GO term: %s<br/> module corrected p-value: %.02e<br/>Module Genes: %s<br/>", go, go_correct_pvalue, paste(unlist(sapply(go_genes, function(x) 
      if(x %in% genename_map$feature_name){
        x <- genename_map$common_name[which(x == genename_map$feature_name)]
      }else{
        x <- x
      })), collapse = ' | ')), tibble())) %>%
    select(module, GO) %>% group_by(module) %>% nest(GO = 'GO') %>% rowwise() %>%
    mutate(GO = paste(unlist(as.list(GO)), collapse = "<br/><br/>"))
    
  
  regulators <- regulators %>% rowwise %>% 
    mutate(Regulators = ifelse(length(module) > 0, sprintf("enriched regulator: %s<br/> regulator corrected p-value: %.02e<br/>regulator targets: %s <br/>", sapply(regulator, function(x)
      if(x %in% genename_map$feature_name){
        x <- genename_map$common_name[which(x == genename_map$feature_name)]
      }else{
        x <- x
      }), reg_correct_pvalue, paste(unlist(sapply(reg_target_genes, function(x) 
      if(x %in% genename_map$feature_name){
        x <- genename_map$common_name[which(x == genename_map$feature_name)]
      }else{
        x <- x
      })), collapse = ' | ')), tibble())) %>%
    select(module, Regulators) %>% nest(Regulators = 'Regulators') %>% rowwise() %>%
    mutate(Regulators = paste(unlist(as.list(Regulators)), collapse = "<br/><br/>"))
  if(method == "list"){
    Module_Table <- Module_Table %>% select(module, "Genes on List", gene_list, corrected_pval) %>% rowwise %>%
      mutate(Genes = paste(sort(unlist(sapply(gene_list, function(x)
        if(x %in% genename_map$feature_name){
          x <- genename_map$common_name[which(x == genename_map$feature_name)]
        }else{
          x <- x
        }))), collapse = ' | ')) %>% 
      mutate("Genes on List" = paste(sort(unlist(sapply(`Genes on List`, function(x)
        if(x %in% genename_map$feature_name){
          x <- genename_map$common_name[which(x == genename_map$feature_name)]
        }else{
          x <- x}))), collapse = ' | ')) %>%
      select(module, "Genes on List", Genes, corrected_pval) %>%
      left_join(GO) %>% left_join(regulators) %>%
      rename("Gene List enrichment p-value" = corrected_pval) %>%
      mutate(module = sprintf('<a href="#" onclick=Shiny.setInputValue("module_id_info", %s);">%s</a>', module, module))
  }else{
  Module_Table <- Module_Table %>% select(module, gene_list) %>% rowwise %>%
    mutate(Genes = paste(sort(unlist(gene_list)), collapse = ' | ')) %>% select(module, Genes) %>%
    left_join(GO) %>% left_join(regulators) %>%
    mutate(module = sprintf('<a href="#" onclick=Shiny.setInputValue("module_id_info", %s);">%s</a>', module, module))
  }
  return(Module_Table)
}


if(file.exists('net_data.Rdata')){
  load('net_data.Rdata')
}



