library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(gridExtra)
library(plotly)


makeSubgraphHeatmapDynamic <- function(subNet,
                                display_name = 1, 
                                #edge_color_by = NA, 
                                #edge_color_palette = "RdBu", 
                                font_size = 18, 
                                direction = 1, 
                                expression_color_palette = "RdBu", 
                                tfa_color_palette = "PiYG",
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
  
  #if(edge_color_by == "Correlation"){
  #  scale_edge_color <- c(-1, 1)
  #}
  
  sym_node_name_by <- ifelse(is.na(node_name_by), NA, sym(node_name_by))
  #sym_edge_color_by <- ifelse(is.na(edge_color_by), NA, sym(edge_color_by))
  
  
  subNet_nodes <- subNet %N>% mutate(component = group_components())  %>% as_tibble() %>% 
    mutate(node_type = ifelse(str_detect(feature, "_nca"), "nca", 
                              ifelse(regulator == "tar", "tar", "scr" ))) 
  
  ### start with NCA
  subset <- subNet_nodes %>% 
    filter(node_type == "nca") %>% 
    mutate(feature = factor(feature, levels = feature)) %>% 
    mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
    unnest_longer(expression)
  
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
        legend.position = "none",
      )
    ggNCA <- ggplotly(ggNCA, tooltip = c("x", "y", "fill"))
  }
  
  
  ## Regulators 
  subset <- subNet_nodes %>% 
    filter(node_type == "scr") %>% 
    mutate(feature = factor(feature, levels = feature)) %>% 
    mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
    unnest_longer(expression)
  
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
    ggReg <- ggplotly(ggReg, tooltip = c("x", "y", "fill"))
  } 
  
  ##Tar
  subset <- subNet_nodes %>% 
    filter(node_type == "tar") %>%
    mutate(feature = factor(feature, levels = feature)) %>% 
    mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
    unnest_longer(expression)
  
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
        plot.margin =margin(t = 0, r = 0, b = 0 , l = 0, unit = 'npc'), 
        legend.position = "none"
      )
    ggTar <- ggplotly(ggTar, tooltip = c("x", "y", "fill"))
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
  
  #gg <- plot_list[[1]]
  #for(i in 2:length(plot_list)){
  #  gg <- gg + plot_list[[i]]
  #}
  
  num_gene_list <- unlist(num_genes_list)
  
  heights <- num_gene_list/sum(num_gene_list)
  gg <- plotly::subplot(
    plot_list,
    nrows = length(heights), 
    heights = heights, 
    shareX = TRUE,
    shareY = TRUE, 
    titleY = TRUE, 
    margin = 0.005
  ) %>% 
  layout(
    showlegend = FALSE,
    coloraxis = list(
      colorbar = list(visible = FALSE)  # <-- hides the colorbar
    )
  )
  
  #row_coords <- lapply(seq_along(num_genes_list), function(i){ 
  #  y_domain <- gg$x$layout[[paste0("yaxis", ifelse(i==1,"",i))]]$domain
  #  n_rows <- num_genes_list[[i]]
    
    # Compute the edges of each row
  #  edges <- seq(y_domain[1], y_domain[2], length.out = n_rows + 1)
    
    # Midpoint of each row = average of top and bottom edge
  #  midpoints <- (edges[-1] + edges[-length(edges)]) / 2
  #  rev(midpoints)
  #})
  
  #row_coords <- unlist(row_coords)
  
  
  #gene_order <- fct_c(fct_rev(nca_gene_order), fct_rev(reg_gene_order), fct_rev(tar_gene_order))
  
  ## Set up network visualization 
  #subNet <- subNet  %N>% mutate(feature = factor(feature, levels(gene_order))) %>% 
  #  arrange(feature) %>%  
  #  mutate(`Common Name` = factor(`Common Name`, levels = `Common Name`)) %>% 
  #  mutate(y = row_coords) %E>%
  #  mutate(x = 0, xend = 0, y = .N()$y[to], yend = .N()$y[from])

  #if(node_name_by == "Common Name"){
  #  subNet <- subNet %E>% 
  #    mutate(edge_name = sprintf('%s -> %s', .N()$`Common Name`[from], .N()$`Common Name`[to]))
  #}else {
  #  subNet <- subNet %E>% 
  #    mutate(edge_name = sprintf('%s -> %s', .N()$feature[from], .N()$feature[to]))
  #}
  
  
  #edges <- subNet %E>% 
  # as_tibble()
  
  #arc_list <- lapply(seq_len(nrow(edges)), function(i){
  #  arc <- create_arc(edges$x[i], edges$y[i], edges$xend[i], edges$yend[i])
  #  arc$edge_name <- edges$edge_name[i]
  #  arc$color_value <- edges[[sym_edge_color_by]][i]
  #  return(arc)
  #})
  #arc_df <- do.call(rbind, arc_list)
  
  #zmin <- min(edges[[sym_edge_color_by]])
  #zmax <- max(edges[[sym_edge_color_by]])
  
  #p <- plot_ly()
  #for(i in seq_len(nrow(edges))) {
  #  arc <- create_arc(edges$x[i], edges$y[i], edges$xend[i], edges$yend[i])
    
    # Map value to a color in a palette
  #  color <- scales::col_numeric(edge_color_palette, domain = c(zmin, zmax), reverse = direction == 1)(edges[[sym_edge_color_by]][i])
    
  #  if (edge_color_by == "Correlation"){
  #    l <- "Correlation"
  #  }else{
  #    l <- "Reg Weight"
  #  }
  #  hover_text <- paste0(edges$edge_name[i], "<br>", l, ": ", sprintf('%.03f', (edges[[sym_edge_color_by]][i])))
    
  #  p <- add_trace(
  #    p,
  #    x = arc$x,
  #    y = arc$y,
  #    type = "scatter",
  #    mode = "lines",
  #    line = list(color = color, width = 2),
  #    text = hover_text,
  #    hoverinfo = "text",
  #    showlegend = FALSE
  #  )
  #}
  
  #p <- p %>% layout(
  #    xaxis = list(
  #    range = c(0, 0.501),
  #    showline = FALSE,
  #    showticklabels = FALSE,
  #    showgrid = FALSE,
  #    zeroline = FALSE
  #  ),
  #  yaxis = list(
  #    range = c(0, 1),
  #    autorange = FALSE,
  #    showline = FALSE,
  #    showticklabels = FALSE,
  #    showgrid = FALSE,
  #    zeroline = FALSE
  #  )
  #)
  
  #gg_combine <- plotly::subplot(gg,
  #                p, 
  #                nrows =1,
  #                widths = c(0.8, 0.2),
  #                margin = 0.
  #)
  
  #gg_combine %>% layout(
  #  plot_bgcolor = "white",
  #  paper_bgcolor = "white"
  #)
  
  return(gg)
}

create_arc <- function(x0, y0, x1, y1, x_offset = -0.5, n = 50) {
  # t goes from 0 to 1
  t <- seq(0, 1, length.out = n)
  
  # y goes linearly from y0 to y1
  y <- y0 + t*(y1 - y0)
  
  x_offset  = abs((y1 - y0) / 2)  
  
  # x is a quadratic that peaks at t=0.5 (middle of the arc)
  # t*(1-t) gives a parabola that is 0 at t=0 and t=1, max at t=0.5
  x <- x0 + (x1 - x0)*t + x_offset*4*t*(1-t)
  
  data.frame(x, y)
}
