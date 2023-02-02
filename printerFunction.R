library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)

makeSubNetGraph <- function(subNet,
                            names_in_nodes = FALSE, node_color_by = NA, 
                            edge_color_by = NA, node_color_palette = NA, 
                            node_size_by = NA, max_node_size = NA, 
                            layout = 'dh', focus_nodes = list(), 
                            font_size = 18, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                            expand_x = 0, expand_y = 0){
  set.seed(88)
  if(layout == "focus"){
    gg <- ggraph(subNet, layout = layout, focus = feature %in% focus_nodes, ) 
  } else{
    gg <- ggraph(subNet, layout = layout) 
  }
  palettes <- tibble(rownames_to_column(brewer.pal.info, var = 'pal'))
 
  ### Edges 
  gg <- gg + 
    theme_light() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
    theme(panel.background = element_rect(fill = "transparent",colour = NA)) + 
    theme(axis.ticks = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.title = element_blank()) 
  
  if(!is.na(edge_color_by)){
    gg <- gg +
      geom_edge_link0(aes_string(color  = edge_color_by))
  }else{
    gg <- gg + 
      geom_edge_link0(color = '#BEBEBE')
  }
  
  
  #### Nodes 
  if(!is.na(node_color_by) && !is.na(node_size_by)){
    gg <- gg + 
      geom_node_point(aes_string(fill = node_color_by, size = node_size_by), shape = 21,  show.legend = show_legend)
  }else if(!is.na(node_color_by)){
    gg <- gg + 
      geom_node_point(aes_string(fill = node_color_by), size = max_node_size, shape = 21, show.legend = show_legend)
  }else if(!is.na(node_size_by)){
    gg <- gg + 
      geom_node_point(aes_string(size = node_size_by), size = max_node_size, fill = 'black', shape = 21, show.legend = show_legend)
  }else{
    gg <- gg +
      geom_node_point()
  }
  
  
  #### Node Labels
  if(names_in_nodes){
    if(!is.na(node_color_by) && !is.na(node_size_by)){
      gg <- gg + 
        geom_node_label(aes_string(fill = node_color_by, size = node_size_by, label = "display_name"), show.legend = FALSE)
    }else if(!is.na(node_color_by)){
      gg <- gg + 
        geom_node_label(aes_string(fill = node_color_by, label = "display_name"), size = font_size, show.legend = FALSE)
    }else if(!is.na(node_size_by)){
      gg <- gg + 
        geom_node_label(aes_string(size = node_size_by, label = "display_name"), show.legend = FALSE)
    }else{
      gg <- gg +
        geom_node_label(aes(label = display_name), size = max_node_size, size = font_size, show.legend = FALSE)
    }
  }else{
    gg<- gg + 
      geom_node_text(aes(label = display_name), size = font_size, nudge_y = nudge_y, angle = text_angle, show.legend = FALSE) 
  }
  
  ### Color Scale
  if(node_color_by == "module"){
    modules <- subNet %N>% pull(module)
    node_breaks <- c("Unlabeled", setdiff(modules, "Unlabeled"))
    palette_val <- palettes %>% filter(pal == node_color_palette ) %>% pull(maxcolors)
    if(length(node_breaks) <= palette_val + 1){
      node_value <- c('#BEBEBE', brewer.pal(max(3, min(length(node_breaks) - 1, palette_val)), node_color_palette))
    }else{
      max_pallette <- brewer.pal(palette_val, node_color_palette)
      extend_pallette <- colorRampPalette(max_pallette)(length(node_breaks) - 1)
      node_value <- c('#BEBEBE', extend_pallette)
    }
  }else if(node_color_by == "geneSuper"){
    geneSuper <- subNet %N>% pull(geneSuper)
    node_breaks <- c('Unlabeled', setdiff(geneSuper, 'Unlabeled'))
    palette_val <- palettes %>% filter(pal == node_color_palette ) %>% pull(maxcolors)
    if(length(node_breaks) <= palette_val + 1){
      node_value <- c('#BEBEBE', brewer.pal(max(3, min(length(node_breaks) - 1, palette_val)), node_color_palette))
    }else{
      max_pallette <- brewer.pal(palette_val, node_color_palette)
      extend_pallette <- colorRampPalette(max_pallette)(length(node_breaks) - 1)
      node_value <- c('#BEBEBE', extend_pallette)
    }
  }else{
    node_breaks <- c('scr', 'tar')
    node_value <- c('#fb8072', '#80b1d3')
  }
  gg <- gg + scale_fill_manual(breaks = node_breaks, 
                     values = node_value)
  
  if(!is.na(edge_color_by)){
    edge_breaks = c(TRUE, FALSE)
    edge_value <- c('#F9B6AF', '#BEBEBE')
    gg <- gg + scale_edge_color_manual(breaks = edge_breaks,
                                       values = edge_value)
  }
  
  if(!is.na(node_size_by)){
    gg <- gg + scale_size_continuous(range = c(1,max_node_size)) 
  }
  
  x_left <- min(gg$data$x)
  x_right <- max(gg$data$x)
  y_bottom <- min(gg$data$y)
  y_top <- max(gg$data$y)
  gg <- gg + coord_cartesian(xlim = c(x_left - expand_x, x_right + expand_x),
                             ylim = c(y_bottom - expand_y, y_top + expand_y))
                             
  return(gg)
}
  
  