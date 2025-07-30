library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(gridExtra)


makeSubgraphHeatmap <- function(subNet,
  display_name = 1, 
  edge_color_by = NA, 
  edge_color_palette = "RdBl", 
  font_size = 18, 
  direction = 1, 
  expression_color_palette = "Reds", 
  tfa_color_palette = "Greens",
  scale_edge_color  = c(-1, 1), 
  scale_expression_colors = c(-5,5), 
  scale_tfa_colors = c(-5, 5), 
  figure_font_size = 18
) 
{ 
  
  if(display_name == 1){
    node_name_by <- "Common Name"
  }else if(display_name == 2){
    node_name_by <- "feature"
  }else{
    print("invalid name type")
    return(ggplot())
  }
  
  
  sym_node_name_by <- ifelse(is.na(node_name_by), NA, sym(node_name_by))
  sym_edge_color_by <- ifelse(is.na(edge_color_by), NA, sym(edge_color_by))
  
  
  subNet_nodes <- subNet %N>% mutate(component = group_components())  %>% as_tibble() %>% 
    mutate(node_type = ifelse(str_detect(feature, "_nca"), "nca", 
                              ifelse(regulator == "tar", "tar", "scr" ))) 
  
  ### start with NCA
  subset <- subNet_nodes %>% 
    filter(node_type == "nca") %>% 
    mutate(feature = factor(feature, levels = feature)) %>% 
    mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
    unnest_longer(expression, indices_to = "expression_id") 
  
  num_nca_genes <- length(unique(subset$feature))
  nca_gene_order <- rev(unique(subset %>% pull(feature)))  
  
  
  if(nrow(subset) > 0){
  
  ggNCA <- ggplot(subset, 
               aes(x = expression_id, y = !!sym_node_name_by, fill = expression)) + 
    geom_tile() + 
    scale_fill_distiller(palette = tfa_color_palette, direction = direction, 
                         limits = scale_tfa_colors, oob = scales::squish, name = "TFA") + 
    scale_x_discrete(expand = c(0,0)) + 
    scale_y_discrete(expand = c(0,0)) + 
    ylab( "TFA") + 
    theme(
      axis.title = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = figure_font_size), 
      legend.text = element_text(size = figure_font_size - 3),
      legend.title = element_text(size = figure_font_size),
      plot.margin = margin(t = 0, r = 0, b = 0 , l = 0, unit = 'npc'),
    ) 
  }
  
  
  ## Regulators 
  subset <- subNet_nodes %>% 
    filter(node_type == "scr") %>% 
    mutate(feature = factor(feature, levels = feature)) %>% 
    mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
    unnest_longer(expression, indices_to = "expression_id") 
  
    num_reg_genes <- length(unique(subset$feature))
    reg_gene_order <- rev(unique(subset %>% pull(feature)))
  
  
  if(nrow(subset) > 0){

    ggReg <- ggplot(subset, 
                    aes(x = expression_id, y = !!sym_node_name_by, fill = expression)) + 
      geom_tile() + 
      scale_fill_distiller(palette = expression_color_palette, direction = direction, 
                           limits = scale_expression_colors, oob = scales::squish) + 
      scale_x_discrete(expand = c(0,0)) + 
      scale_y_discrete(expand = c(0,0)) + 
      ylab("Regulators") + 
      theme(
        axis.title = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = figure_font_size), 
        legend.text = element_text(size = figure_font_size - 3),
        legend.title = element_text(size = figure_font_size),
        plot.margin = margin(t = 0, r = 0, b = 0 , l = 0, unit = 'npc'),
        legend.position = "none"
      )
  } 
  
  ##Tar
  subset <- subNet_nodes %>% 
    filter(node_type == "tar") %>%
    mutate(feature = factor(feature, levels = feature)) %>% 
    mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
    unnest_longer(expression, indices_to = "expression_id") 
  
  num_tar_genes <- length(unique(subset$feature))
  tar_gene_order <- rev(unique(subset %>% pull(feature)))
 
  if(nrow(subset) > 0){
  
    ggTar <- ggplot(subset, 
                    aes(x = expression_id, y = !!sym_node_name_by, fill = expression)) + 
      geom_tile() + 
      scale_fill_distiller(palette = expression_color_palette, direction = direction, 
                           limits = scale_expression_colors, oob = scales::squish) + 
      scale_x_discrete(expand = c(0,0)) + 
      scale_y_discrete(expand = c(0,0)) + 
      ylab("Targets") + 
      theme(
        axis.title = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = figure_font_size), 
        legend.text = element_text(size = figure_font_size - 3),
        legend.title = element_text(size = figure_font_size),
        plot.margin =margin(t = 0, r = 0, b = 0 , l = 0, unit = 'npc')
      )
  }
  
  ## Make combine plot_list
  plot_list <- list()
  num_genes_list <- list()
  
  if(num_nca_genes > 0){
    plot_list  <- append(plot_list, list(ggNCA))
    num_genes_list <- append(num_genes_list, num_nca_genes)
  }
  
  if(num_reg_genes > 0){
    plot_list <- append(plot_list, list(ggReg))
    num_genes_list <- append(num_genes_list, num_reg_genes)
  }
  
  if(num_tar_genes > 0){
    plot_list <- append(plot_list, list(ggTar))
    num_genes_list <- append(num_genes_list, num_tar_genes)
  }
  
  gg <- plot_list[[1]]
  for(i in 2:length(plot_list)){
    gg <- gg + plot_list[[i]]
  }
  
  gg <- gg + plot_layout(ncol = 1, nrow = length(plot_list), heights = num_genes_list,guides = 'collect')

  gene_order <- fct_c(tar_gene_order, reg_gene_order, tar_gene_order)
    
  ## Set up network visualization 
  subNet <- subNet  %N>% mutate(feature = factor(feature, levels(gene_order))) %>% 
    arrange(feature) %>%  
    mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
    mutate(y = 0.5:1:nrow(subNet %N>% as_tibble())) #%E>% 
    #mutate(forward_alpha = ifelse(.N()$y[from] > .N()$y[to], 0, 1)) %>% 
    #mutate(reverse_alpha = ifelse(.N()$y[from] < .N()$y[to], 0, 1))
  
  
  net <- ggraph(subNet, x = 0, y = y, layout = "manual") +
    geom_edge_arc(aes(color = !!sym_edge_color_by), arrow = arrow(angle = 15, ends ='last', length = unit(0.15, "inches"), type = 'closed'), strength = -0.05) + 
    geom_edge_arc(aes(color = !!sym_edge_color_by), arrow = arrow(angle = 15, ends ='last', length = unit(0.15, "inches"), type = 'closed'), strength = 0.05) +  
    #geom_node_point() + 
    scale_y_continuous(expand = c(0,0), limits = c(0, max(subNet %N>% as_tibble %>% pull(y)) + .5)) + 
    scale_x_continuous(expand = c(.05,.05), limits = c(NA, 0)) + 
  
    scale_edge_color_distiller(palette = edge_color_palette, direction = -1, limits = scale_edge_color, oob = scales::squish) + 
    scale_alpha_continuous(limits = c(0, 1), range= c(0, 1)) + 
    theme_void() + 
    theme(
      legend.position = "left",
      legend.text = element_text(size = figure_font_size - 3),
      legend.title = element_text(size = figure_font_size),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'npc')
    )
    
  gg_combine <- (net | gg) + 
    plot_layout(ncol = 2, widths = c(1, 5))

  return(gg_combine)
}
