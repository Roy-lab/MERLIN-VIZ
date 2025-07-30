library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(ggpmisc)

makeSubNetGraph <- function(subNet,
                            names_in_nodes = FALSE, node_color_by = NA, 
                            edge_color_by = NA, edge_color_palette = "RdBu", 
                            node_color_palette = 'Dark2', 
                            node_size_by = NA, max_node_size = NA,
                            edge_width = 1, 
                            layout = 'dh', focus_nodes = list(), 
                            font_size = 18, nudge_y = 0, text_angle  = 0, show_legend = TRUE, direction = 1,
                            expand_x = 0, expand_y = 0, font_color = '#ffffff', unlab_color = '#000000',
                            arrow_size = 0.5, 
                            color_scale_limits = c(-5,5) , legend_font_size = 18)
  {
  set.seed(88)
  sym_node_color_by <- ifelse(is.na(node_color_by), NA, sym(node_color_by))
  sym_node_size_by <-  ifelse(is.na(node_size_by), NA, sym(node_size_by))
  sym_edge_color_by <- ifelse(is.na(edge_color_by), NA, sym(edge_color_by))
  #print(edge_color_by)
  subNet <- subNet %>% mutate(display_name = ifelse(is.na(display_name), yes = "", no = display_name))
  
  if(layout == "focus"){
    gg <- ggraph(subNet, layout = layout, focus = feature %in% focus_nodes) 
  } else{
    gg <- ggraph(subNet, layout = layout) 
  }
  palettes <- tibble(rownames_to_column(brewer.pal.info, var = 'pal'))
  
  
  gg <- gg + 
    theme_light() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
    theme(panel.background = element_rect(fill = "transparent",colour = NA)) + 
    theme(axis.ticks = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.title = element_blank()) 
  
  
  
  ### Edges
  if(!is.na(edge_color_by)){
    gg <- gg +
      geom_edge_link(aes(color  = !!sym_edge_color_by, 
                         start_cap = label_rect(node1.display_name, fontsize = font_size), 
                         end_cap = label_rect(node2.display_name, fontsize =font_size)), edge_width = edge_width, 
                     arrow = arrow(angle = 15, ends ='last', length = unit(arrow_size,  "lines"), type = 'closed'), 
                     #end_cap =  circle(2, 'mm'),
                     show.legend = TRUE)
  }else{
    gg <- gg + 
      geom_edge_link(aes(start_cap = label_rect(node1.display_name, fontsize = font_size), 
                         end_cap = label_rect(node2.display_name, fontsize = font_size)),
                     color = '#BEBEBE', 
                     arrow = arrow(ends ='last', angle = 15, length = unit(arrow_size, "lines"), type = 'closed'),
                     #end_cap = rectangle(width = 10, height = 10, width_unit = 'mm', height_unit = 'mm'),
                     edge_width = edge_width)
  }
  
  
  #### Nodes 
  if(!is.na(node_color_by) && !is.na(node_size_by)){
    gg <- gg + 
      geom_node_point(aes(fill = !!sym_node_color_by, shape = regulator, size = !!sym_node_size_by),  show.legend = FALSE)
  }else if(!is.na(node_color_by)){
    gg <- gg + 
      geom_node_point(aes(fill = !!sym_node_color_by, shape = regulator), size = max_node_size, show.legend = FALSE)
  }else if(!is.na(node_size_by)){
    gg <- gg + 
      geom_node_point(aes(size = !!sym_node_size_by, shape = regulator), size = max_node_size, fill = 'black', show.legend = FALSE)
  }else{
    gg <- gg +
      geom_node_point()
    
  }
  
  
  #### Node Labels
  if(names_in_nodes){
    if(!is.na(node_color_by) && !is.na(node_size_by)){
      gg <- gg + 
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), fill = !!sym_node_color_by, size = !!sym_node_size_by, label = display_name), color = font_color, label.r = unit(0, "pt"), show.legend = FALSE) + 
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), fill = !!sym_node_color_by, label = display_name), label.r = unit(.25, "lines"), color = font_color, size = font_size, show.legend = FALSE)
    }else if(!is.na(node_color_by)){
      gg <- gg + 
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), fill = !!sym_node_color_by, label = display_name), color = font_color, label.r = unit(0, "lines"), size = font_size, show.legend = FALSE) +
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), fill = !!sym_node_color_by, label = display_name), color = font_color, label.r = unit(.25, "lines"), size = font_size, show.legend = FALSE)
    }else if(!is.na(node_size_by)){
      gg <- gg + 
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), size = !!sym_node_size_by, label = display_name), color = font_color, label.r = unit(0, "lines"), show.legend = FALSE) + 
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), size = !!sym_node_size_by, label = display_name), color = font_color, label.r = unit(.25, "lines"), show.legend = FALSE)
    }else{
      gg <- gg +
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), label = display_name), label.r = unit(0, "lines"), color = font_color, size = max_node_size, size = font_size, show.legend = FALSE) + 
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), label = display_name), label.r = unit(0, "lines"), color = font_color, size = max_node_size, size = font_size, show.legend = FALSE)
    }
  }else{
    gg<- gg + 
      geom_node_text(aes(filter = display_name != "", label = display_name), size = font_size, nudge_y = nudge_y, angle = text_angle, show.legend = FALSE) 
  }
  
  
  
  ### Color Scale
  if(node_color_by == "module"){
    modules <- subNet %N>% pull(module)
    node_breaks <- c("Unlabeled", setdiff(modules, "Unlabeled"))
    palette_val <- palettes %>% filter(pal == node_color_palette ) %>% pull(maxcolors)
    if(length(node_breaks) <= palette_val + 1){
      node_value <- c(unlab_color, brewer.pal(max(3, min(length(node_breaks) - 1, palette_val)), node_color_palette))
    }else{
      max_palette <- brewer.pal(palette_val, node_color_palette)
      extend_palette <- colorRampPalette(max_palette)(length(node_breaks))
      node_value <- c(unlab_color, extend_palette)
    }
  }else if(node_color_by == "geneSuper"){
    geneSuper <- subNet %N>% pull(geneSuper)
    node_breaks <- c('Unlabeled', setdiff(geneSuper, 'Unlabeled'))
    palette_val <- palettes %>% filter(pal == node_color_palette ) %>% pull(maxcolors)
    if(length(node_breaks) <= palette_val + 1){
      node_value <- c(unlab_color, brewer.pal(max(3, min(length(node_breaks) - 1, palette_val)), node_color_palette))
    }else{
      max_palette <- brewer.pal(palette_val, node_color_palette)
      extend_palette <- colorRampPalette(max_palette)(length(node_breaks))
      node_value <- c(unlab_color, extend_palette)
    }
  }else{
    node_breaks <- c('scr', 'tar')
    node_value <- c('#fb8072', '#80b1d3')
  }
  gg <- gg + scale_fill_manual(breaks = node_breaks, values = node_value)
  
  
  
  if(edge_color_by == "is_steiner"){
    edge_breaks = c(TRUE, FALSE)
    edge_value <- c('#fb8072', '#BEBEBE')
    gg <- gg + scale_edge_color_manual(breaks = edge_breaks, values = edge_value)
  }else if(edge_color_by == "Correlation"){
    gg <- gg + scale_edge_color_distiller(palette = edge_color_palette, direction = -1, limits = color_scale_limits, oob = scales::squish)
  }else if(edge_color_by == "Reg_weight") {
    gg <- gg + scale_edge_color_distiller(palette = edge_color_palette, direction = -1, limits = color_scale_limits, oob = scales::squish)
  }
  
  if(!is.na(node_size_by)){
    gg <- gg + scale_size_continuous(range = c(1,max_node_size)) 
  }
  gg <- gg + scale_shape_manual(breaks = c('scr', 'tar'), values =c(22, 21)) + 
    guides(fill=guide_legend(override.aes=list(shape=21)))
  
  x_left <- min(gg$data$x)
  x_right <- max(gg$data$x)
  y_bottom <- min(gg$data$y)
  y_top <- max(gg$data$y) 
  gg <- gg + coord_cartesian(xlim = c(x_left - expand_x, x_right + expand_x),
                             ylim = c(y_bottom - expand_y, y_top + expand_y)) 
  
  if(edge_color_by == "Correlation")
  {
    title = "Correlation"
  } else if ( edge_color_by == "Reg_weight")
  {
    title = "Regression weight"
  } else {
    title = "Steiner"
  }
  
  gg <- gg + 
    theme(legend.position  =  "bottom", 
          legend.title = element_text(size = legend_font_size), 
          legend.text = element_text(size = legend_font_size - 3, color = "black"), 
          legend.background = element_rect(fill = "white")) + 
    guides(
      #\edge_color = guide_colorbar(title = title),  # Keep the colorbar legend for edges
      edge_linetype = "none",                     # You should handle this in the geom
      edge_arrow = "none",                        # Arrows should also be handled in the geom
      shape = "none",                             # Remove shape legend
      fill = "none",                              # Remove fill legend
      color = "none",                            # Remove other color legends
    )
    
  
  return(gg)
}

