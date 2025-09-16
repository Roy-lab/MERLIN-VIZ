library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(ggpmisc)

makeSubNetGraph <- function(subNet,
                            names_in_nodes = FALSE, node_color_by = NA, 
                            edge_color_by = NA,  
                            edge_width_by = NA, 
                            edge_color_palette = "RdBu", 
                            node_color_palette = 'Dark2', 
                            max_edge_width = NA, 
                            node_size_by = NA, max_node_size = 5,
                            edge_width = 1, 
                            layout = 'dh', focus_nodes = list(), 
                            font_size = 18, nudge_y = 0, text_angle  = 0, show_legend = TRUE, direction = 1,
                            expand_x = 0, expand_y = 0, font_color = '#ffffff', unlab_color = '#000000', 
                            exp_scale_limits = c(0, 5), 
                            color_scale_limits = c(-5,5) , legend_font_size = 18)
  {
  
  ## initial symbol objects ----
  set.seed(88)
  sym_node_color_by <- ifelse(is.na(node_color_by), NA, sym(node_color_by))
  sym_node_size_by <-  ifelse(is.na(node_size_by), NA, sym(node_size_by))
  
  sym_edge_color_by <-sym(edge_color_by)
  sym_edge_width_by <- sym(edge_width_by)
  
  ## fill NA display names with empty -----
  subNet <- subNet %>% mutate(display_name = ifelse(is.na(display_name), yes = "", no = display_name))
  
  ## Initialize layout ----
  if(layout == "focus"){
    gg <- ggraph(subNet, layout = layout, focus = feature %in% focus_nodes) 
  } else{
    gg <- ggraph(subNet, layout = layout) 
  }
  
  ## Initialize palettes ---- 
  palettes <- tibble(rownames_to_column(brewer.pal.info, var = 'pal'))
  
  
  ## Set ggplot background ----
  gg <- gg + 
    theme_light() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + 
    theme(panel.background = element_rect(fill = "transparent",colour = NA)) + 
    theme(axis.ticks = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.title = element_blank()) 
  
  
  
  ## Edges display ---- 
  gg <- gg +
      geom_edge_link(aes(color  = !!sym_edge_color_by, 
                         start_cap = label_rect(node1.display_name, fontsize = font_size), 
                         end_cap = label_rect(node2.display_name, fontsize =font_size), edge_width = abs(!!sym_edge_width_by)), 
                     arrow = arrow(angle = 15, ends ='last', length = unit(.5, "lines"), type = 'closed'), 
                     #end_cap =  circle(2, 'mm'),
                     show.legend = TRUE)
  
  
  ## Nodes display ---- 
  ### Node with color and size ----
  if(!is.na(node_color_by) && !is.na(node_size_by)){
    gg <- gg + 
      geom_node_point(aes(fill = !!sym_node_color_by, shape = regulator, size = !!sym_node_size_by)) #,  show.legend = FALSE)
  ## node with color, default size ---- 
  }else if(!is.na(node_color_by)){
    gg <- gg + 
      geom_node_point(aes(fill = !!sym_node_color_by, shape = regulator), size = max_node_size) # , show.legend = FALSE)
  ## node with size, default color (black) -----
  }else if(!is.na(node_size_by)){
    gg <- gg + 
      geom_node_point(aes(size = !!sym_node_size_by, shape = regulator), size = max_node_size, fill = 'black') #, show.legend = FALSE)
  ## Default node -----  
  }else{
    gg <- gg +
      geom_node_point()
  }
  
  
  ## Node Labels ----
  ### Node label within plot ----
  if(names_in_nodes){
    #### Node with color and size -----
    if(!is.na(node_color_by) && !is.na(node_size_by)){
      gg <- gg + 
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), fill = !!sym_node_color_by, size = !!sym_node_size_by, label = display_name), color = font_color, label.r = unit(0, "pt"), show.legend = FALSE) + 
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), fill = !!sym_node_color_by, label = display_name), label.r = unit(.25, "lines"), color = font_color, size = font_size, show.legend = FALSE)
    #### Node with color ----  
    }else if(!is.na(node_color_by)){
      gg <- gg + 
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), fill = !!sym_node_color_by, label = display_name), color = font_color, label.r = unit(0, "lines"), size = font_size, show.legend = FALSE) +
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), fill = !!sym_node_color_by, label = display_name), color = font_color, label.r = unit(.25, "lines"), size = font_size, show.legend = FALSE)
    #### Node with size ----
    }else if(!is.na(node_size_by)){
      gg <- gg + 
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), size = !!sym_node_size_by, label = display_name), color = font_color, label.r = unit(0, "lines"), show.legend = FALSE) + 
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), size = !!sym_node_size_by, label = display_name), color = font_color, label.r = unit(.25, "lines"), show.legend = FALSE)
    #### Node default with labels ----  
    }else{
      gg <- gg +
        geom_node_label(aes(filter = (regulator == 'scr' & display_name != ""), label = display_name), label.r = unit(0, "lines"), color = font_color, size = max_node_size, size = font_size, show.legend = FALSE) + 
        geom_node_label(aes(filter = (regulator == 'tar' & display_name != ""), label = display_name), label.r = unit(0, "lines"), color = font_color, size = max_node_size, size = font_size, show.legend = FALSE)
    }
    
    ## Node name outside (default) ----
  }else{
    if(!is.na(node_size_by)){
      gg<- gg + 
        geom_node_text(aes(filter = display_name != "", label = display_name, size = !!sym_node_size_by), nudge_y = nudge_y, angle = text_angle, show.legend = FALSE) 
    }else{
      gg<- gg + 
        geom_node_text(aes(filter = display_name != "", label = display_name), size = font_size, nudge_y = nudge_y, angle = text_angle, show.legend = FALSE) 
    }
  }
  
  
  
  ## Node color scale ----
  skip_discrete_palette <- FALSE
  ### Module color scale (discrete) -----
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
  ### Gene Super (first three character) color scale (discrete) ----
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
  ### Node by gene expression value (continuous)  ----
  }else if(node_color_by %in% c("exp","mean_expression")){
    gg <- gg + scale_fill_distiller(palette = node_color_palette, oob = scales::squish, direction = 1)
    skip_discrete_palette <- TRUE
  ### Node by node type (discrete) ------  
  }else{
    node_breaks <- c('scr', 'tar')
    node_value <- c('#fb8072', '#80b1d3')
  }
  
  
  #### discrete method scale implementation  -----
  if(!skip_discrete_palette){
    gg <- gg + scale_fill_manual(breaks = node_breaks, values = node_value)
  }
  
  ## Edge color scale  ---- 
  ### Edge steiner -----
  if(edge_color_by == "is_steiner"){
    edge_breaks = c(TRUE, FALSE)
    edge_value <- c('#fb8072', '#BEBEBE')
    gg <- gg + scale_edge_color_manual(breaks = edge_breaks, values = edge_value)
  ### Correlation ----  
  #}else if(edge_color_by == "Correlation"){
  #  gg <- gg + scale_edge_color_distiller(palette = edge_color_palette, direction = -1, limits = c(-1, 1))
  ### Regression ----
  #}else if(edge_color_by == "Reg_weight") {
  #  gg <- gg + scale_edge_color_distiller(palette = edge_color_palette, direction = -1, limits = color_scale_limits, oob = scales::squish)
  }else{
    gg <- gg + scale_edge_color_distiller(palette = edge_color_palette, direction = -1, limits = color_scale_limits, oob = scales::squish, name = "Regression Weight")
  }
  
  ## Node size scale ----
  if(!is.na(node_size_by)){
    gg <- gg + scale_size_continuous(range = c(1,max_node_size)) 
  }
  
  ## Node shape scale ----
  gg <- gg + scale_shape_manual(breaks = c('scr', 'tar'), values =c(23, 22))
  
  ## edge width scale -----
  gg <- gg + scale_edge_width(limits = c(0, 1), range = c(0, max_edge_width), name = 'Absolute Correlation')
  
  
  ## Plot plane limit setup (expand) ----
  x_left <- min(gg$data$x)
  x_right <- max(gg$data$x)
  y_bottom <- min(gg$data$y)
  y_top <- max(gg$data$y) 
  gg <- gg + coord_cartesian(xlim = c(x_left - expand_x, x_right + expand_x),
                             ylim = c(y_bottom - expand_y, y_top + expand_y)) 
  
  ## Title setup for edge ----
  # if(edge_color_by == "Correlation")
  # {
  #   title = "Correlation"
  # } else if ( edge_color_by == "Reg_weight")
  # {
  #   title = "Regression weight"
  # } else {
  #   title = "Steiner"
  # }
  title = "Regression weight"
  
  ## Legend initialization -----
  gg <- gg + 
    theme(legend.position  =  "bottom", 
          legend.title = element_text(size = legend_font_size), 
          legend.text = element_text(size = legend_font_size - 3, color = "black"), 
          legend.background = element_rect(fill = "white")) + 
    ### legend option removal ---- 
    guides(
      #\edge_color = guide_colorbar(title = title),  # Keep the colorbar legend for edges
      #size = "none", 
      edge_linetype = "none",                     # You should handle this in the geom
      edge_arrow = "none",                        # Arrows should also be handled in the geom
      #shape = "none",                             # Remove shape legend
      #fill = "none",                              # Remove fill legend
      #color = "none",                            # Remove other color legends
    )
    
  
  return(gg)
}

