#spencers new version that im changing but now adding new tabs

library(shiny)
library(tools)
library(tidyverse)
library(Matrix)
library(scales)
library(DT)
library(webshot)
library(htmlwidgets)
library(RColorBrewer)
library(shinyWidgets)
library(shinythemes)
library(shinyBS)
library(htmltools)
library(bsplus)
library(shinyjs)
library(patchwork)


## Load in aux functions and data
source('aux_functions.R') 
source('printerFunction.R')
source('heatmap_printer.R')
source('heatmap_printer_dynamic.R')
source('ui_items.R')

orthologs_1_to_1 <- Net %N>% pull(`Ortholog 1-1`) %>% unique() %>% .[!is.na(.)]  

ortholog_map <- Net %N>% 
  as_tibble() %>% 
  select(feature, `Ortholog 1-1`) %>% 
  filter(!is.na(`Ortholog 1-1`))

all_gene_names <- unique(c(genes, genename_map$common_name, orthologs_1_to_1))

module_ids <- setdiff(Net %N>% pull(module), c(-9999))
palettes_nodes<- tibble(rownames_to_column(brewer.pal.info, var = 'pal')) %>% 
  filter(category == "qual")
palettes_edges <- tibble(rownames_to_column(brewer.pal.info, var = 'pal')) %>% 
  filter(category == "div")
igraph_layout <- c('Fruchterman-Reingold'='nicely', 'Davidson-Harel'='dh', 'Kamada-Kawai'='kk', 'Large graph layout'= 'lgl') #'Force directed' = 'drl')

################################# ui ###########################################
ui <- navbarPage("GRAsp",
                 id = 'navbar',
                 theme = shinytheme("flatly"), #CC: Set this as the main theme
                 
                 tabPanel("Visualize",
                          #shinythemes::themeSelector(), #CC: this allows you to cycle through different themes when the app is running
                          #####  Gene selection Options ##############################
                          fluidRow(
                            column(2,
                                   method_Picker(),
                                   #### Module Selection ###################
                                   conditionalPanel(
                                     condition = "input.method == 'module'",
                                     module_Picker(module_ids)
                                   ),
                                   #### GO Selection  ###################
                                   conditionalPanel(
                                     condition = "input.method == 'go_term'",
                                     pickerInput(inputId = "go_term", label = "GO terms", choices = c("", unlist(sort(enriched_go_terms))))
                                   ),
                                   
                                   ### Gene List Selection ###############
                                   conditionalPanel(
                                     condition = "input.method == 'list'",
                                     geneList_Picker_Tag(),
                                     geneList_Picker(),
                                     geneList_File_Tag(),
                                     geneList_File()
                                   ),
                                   
                                   ############# Gene List Selection  2 ###########################
                                   conditionalPanel(
                                     condition = "input.method == 'list'",
                                     geneList_checkBox()
                                   ),
                                   ####### Diffusion Selection ##################
                                   conditionalPanel(
                                     condition = "input.method == 'diff'",
                                     diffusion_File_Tag(), 
                                     diffusion_File(), 
                                     diffusion_minTarget_Numerical(), 
                                     diffusion_Lambda_Select(), 
                                     diffusion_topRegulator_Numeric(), 
                                     diffusion_refresh_action()
                                   ),
                            ),
                            
                            #### Visualization Block ########
                            column(10,
                                   tabsetPanel(id = "displayType", type = "tabs",
                                               tabPanel("Network Plot", plotOutput("print_net", click = 'plot_click', height = '1000px')), #CC: Changed to more descriptive titles
                                               tabPanel("Expression Heatmaps", plotlyOutput("expression_heatmap", height = "1000px")),
                                               tabPanel("Gene Table", DT::dataTableOutput("nodes_table")),
                                               tabPanel("Module Table", DT::dataTableOutput("module_table"))
                                               
                                   )
                            )
                          ),
                          fluidRow(
                            ############# Network Plot support options  ########################
                            conditionalPanel(
                              condition ="input.displayType == 'Network Plot'",
                              column(2),
                              column(2,
                                     networkViz_layout_Select(igraph_layout),
                                     networkViz_minGeneCC_Slider(),
                                     networkViz_dispName_Select(),
                                     networkViz_nameFormat_Radio(),
                                     networkViz_namePositionBool_checkBox(),
                                     conditionalPanel(
                                       condition = "input.print_name_bool == 0",
                                       networkViz_nudgeY_Slider(),
                                       networkViz_nameAngle_Slider(), 
                                     )
                              ),
                              column(2,
                                     networkViz_nodeColor_Radio(),
                                     conditionalPanel(
                                       condition = "input.print_group_by == 'module'",
                                       networkViz_nodeColor_Select(palettes_nodes)
                                     ),
                                     networkViz_nodeSize_Slider(),
                                     networkViz_nodeFontSize_Slider()
                                     
                              ),
                              column(2,
                                     networkViz_edgeColor_Radio(),
                                     conditionalPanel(
                                       condition = "input.edge_color_by == 'Reg_weight'",
                                       networkViz_edgeRangeReg_Slider()
                                     ),
                                     networkViz_edgePalette_Select(palettes_edges)
                              ),
                              column(2,
                                     networkViz_expandX_Slider(), 
                                     networkViz_expandY_Slider(),
                                     networkViz_legendSize_Slider(), 
                                     networkViz_imageSaveHeight_Slider(), 
                                     networkViz_imageSaveWidth_Slider(), 
                                     networkViz_imageName_Text(),
                                     networkViz_fileType_Select(), 
                                     networkViz_Download()
                              )
                            ),
                            ############# Heatmap Support option Tab  ########################
                            conditionalPanel(
                              condition ="input.displayType == 'Expression Heatmaps'",
                              column(2),
                              column(2,
                                     heatmapViz_nameFormat_Radio(),
                              ),
                              column(2,
                                     heatmapViz_TFAPalette_Select(palettes_edges),
                                     heatmapViz_TFARange_Slider(),
                                     heatmapViz_expPalette_Select(palettes_edges),
                                     heatmapViz_expRange_Slider() 
                              ),
                              # column(2,
                              #        heatmapGraphViz_edgeColor_Radio(),
                              #        conditionalPanel(
                              #          condition = "input.edge_color_by_heatmap == 'Reg_weight'",
                              #          heatmapGraphViz_edgeRangeReg_Slider(),
                              #        ),
                              #        heatmapGraphViz_edgePalette_Select()
                              # ),
                              column(2,
                                     heatmapViz_fontSize_Slider(), 
                                     heatmapViz_imageHeight_Slider(),
                                     heatmapViz_imageWidth_Slider(),
                                     heatmapViz_fileName_Text(), 
                                     heatmapViz_fileType_Select(),
                                     heatmapViz_Download()
                              )
                            )
                          )
                 ),
                 
                 ######### 
                 
                 tabPanel("About",
                          fluidPage(htmltools::tags$iframe(src = "Help.html", width = '100%', height = 1000, style = "border:none;"))
                 ),
                 tabPanel("Contact", #CC: connected a googleforms that allows users to ask questions or give suggestions. Not sure if there was a better way.
                          tags$iframe(
                            src="https://docs.google.com/forms/d/e/1FAIpQLScwQKwts37i7ykKs1-wlpcE-fPOVVyjUYzLBXTdE93vI6zPLA/viewform?embedded=true",
                            width = "100%",
                            height = "600",
                            frameborder = "0",
                            marginheight = "0",
                            marginwidth = "0",
                            scrolling = "auto"
                        )),
                 
                # JavaScript to switch tabs
                tags$script(HTML("
                Shiny.addCustomMessageHandler('switchTab', function(tabName) {
                var tabLink = $('a:contains(\"' + tabName + '\")');
                if (tabLink.length) tabLink.click();
                });
                "))
                 
                        
)
##################### Server Functions ########################################
server <- function(input, output, session) {
  updateSelectizeInput(session, 'gene', choices = all_gene_names, server = TRUE)
  updateSelectizeInput(session, 'gl', choices = all_gene_names, selected = "srbA", server = TRUE)
 
  
  
  
  
  #### Variables ################
  node_name_info <- reactiveVal(value = NA)  
  module_id_info <- reactiveVal(value = NA)
  steiner_net <- reactiveVal()
  sub_net <- reactiveVal()
  gene_name <-reactiveVal()
  percentile <- reactiveVal(95)
  min_neigh <- reactiveVal(5)
  disp_regs <- reactiveVal(5)
  diff_nodes <- reactiveVal()
  render_diff <- reactiveVal(FALSE)
  gg_out_plot <- reactiveVal(NULL)
  gg_out_heatmap <- reactiveVal(NULL)
  disp_nodes <- reactiveVal(NULL)
  
  lambda <- reactiveVal(1)
  gene_list <- reactiveVal(NULL)
  
  
  ##Contact us action item 
  #bserveEvent(input$contact_button, {
  #  js$windowOpen("https://forms.gle/Xf9S35TatoW3Vmgx9", "_blank")
  #})
  
  ##### File IO ################
  file_path <- reactive({
    file <- input$cell_list_file
    #print(file)
    #req(file)
    file$datapath 
  })
  
  diff_file_path <- reactive({
    file<-input$diff_list_file
    req(file)
    file$datapath
  })
  
  observeEvent(input$refresh_diff, {
    percentile(input$percentile)
    #print(percentile())
    lambda(input$kernel)
    #print(lambda())
  })
  
  observeEvent(input$min_neigh, {
    min_neigh(input$min_neigh)
    render_diff(TRUE)
  })
  
  observeEvent(input$disp_regs, {
    disp_regs(input$disp_regs)
    render_diff(TRUE)
  })
  
  observeEvent(render_diff(), {
    nodes <- diff_nodes()
    if(!is.null(nodes$score)){
      Net <- left_join(Net, nodes)
      sub_net(diffScoreSubgraph(Net, min_neigh(), disp_regs()))
    }
    render_diff(FALSE)
    #print(render_diff())
  })
  
  observeEvent(input$gl,{
    #print(length(Net %N>% pull(feature)))
    #print(input$gl)
    if(length(input$gl) == length(Net %N>% pull(feature)) | length(input$gl) == 0){
      gene_list(NULL)
    }
    else{
      temp <- input$gl
      g <- sapply(temp, function(x) 
        if(x %in% genename_map$common_name){
          x <- genename_map$feature_name[which(x == genename_map$common_name)]
        }
        else if( x %in% ortholog_map$`Ortholog 1-1`){
          x <- ortholog_map$feature[which(x == ortholog_map$`Ortholog 1-1`)]
        }
        else{
          x <- x
        })
      gene_list(g)
    }
  }) 
  
  output$fileUploaded <- reactive({
    return(!is.null(sub_net()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  output$diffFileUploaded <- reactive({
    return(!is.null(sub_net()))
  })
  outputOptions(output, 'diffFileUploaded', suspendWhenHidden=FALSE)
  
  observeEvent(input$go_term,{
    if(input$go_term == ""){
      sub_net(tbl_graph())	
    }else{
      go_term <- input$go_term
      #print(go_term)
      subNet <- goSubgraph(Net, Module, enrich_2_module, go_term)
      Nodes <- subNet %N>% as_tibble()
      regulators <- subNet %N>% as_tibble() %>%  filter(regulator == 'scr') 
      if(input$common_name == 1){
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`),  selected = regulators %>% pull(`Common Name`), server = TRUE)
      }else{
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(feature), selected = regulators %>% pull(feature), server = TRUE)
      }
      sub_net(subNet)
    }
  })
  
  observeEvent(input$module_id,{
    #print(input$module_id)
    if(input$module_id  == ""){
      sub_net(tbl_graph())
    }else{
      module <- input$module_id
      #print(module)
      module_id_info(input$module_id)
      node_name_info(NA)
      subNet <- moduleSubgraph(Net, Module, module)
      Nodes <- subNet %N>% as_tibble()
      regulators <- subNet %N>% as_tibble() %>%  filter(regulator == 'scr') 
      
      
      if(input$common_name == 1){
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`),  selected = regulators %>% pull(`Common Name`), server = TRUE)
      }else{
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(feature), selected = regulators %>% pull(feature), server = TRUE)
      }
      sub_net(subNet)
      
    }
  })
  
  observeEvent(diff_file_path(), {
    id <- showNotification("Computing Defused Scores...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    fp <- diff_file_path()
    score_list <- read_csv(file = fp, col_names = c("feature", "score"))
    
    if(all(is.na(score_list$score))){
      score_list<- score_list %>% mutate(score = 100)
    }
    if(lambda() == 1){
      if(!"k1_sparse" %in% ls()){
        load('k1.Rdata')
      }
      Net <- computeDiffusionScore(Net, score_list, k1_sparse)
    }else if(lambda() ==10){
      if(!"k10_sparse" %in% ls()){
        load('k10.Rdata')
      }
      Net <- computeDiffusionScore(Net, score_list, k10_sparse)
    }else if(lambda() ==100){
      if(!"k100_sparse" %in% ls()){
        load('k100.Rdata')
      }
      Net <- computeDiffusionScore(Net, score_list, k100_sparse)
    }else if(lambda() ==1000){
      if(!"k1000_sparse" %in% ls()){
        load('k1000.Rdata')
      }
      Net <- computeDiffusionScore(Net, score_list, k1000_sparse)
    }
    diff_nodes(Net %N>% as_tibble())
    #print("HERE:::")
    render_diff(TRUE)
  })
  
  observeEvent(file_path(), {
    fp <- file_path()
    sub_net(gene_list(read_csv(file = fp, col_names=FALSE) %>% pull(X1)))
  })
  
  observeEvent(gene_list(), {
    if(length(gene_list() > 0 )){
      subNet <- geneListSubgraph(Net, Module, gene_list(), input$search_additional)
      Nodes <- subNet %N>% as_tibble()
      Nodes_gene_list <- Nodes %>% filter(feature %in% gene_list()) 
      
      if(init){
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`), selected = disp_names, server = TRUE)
        init <<- FALSE
      }else if(input$common_name == 1){
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`),  selected = Nodes_gene_list %>% pull(`Common Name`), server = TRUE)
      }else{
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(feature), selected = Nodes_gene_list %>% pull(feature), server = TRUE)
      }
      sub_net(subNet)
    }else{
      updateSelectizeInput(session, 'print_disp_names', choices = c(""),  selected = c(""), server = TRUE)
    }
  })
  
  
  observeEvent(input$common_name, {
    if(init == TRUE){
      return()
    }
    SB <- sub_net() 
    nodes <- SB %N>% as_tibble()
    
    if(nrow(nodes) == 0)
    {
      updateSelectizeInput(session, 'print_disp_names', choices = c(""), selected = c(""), server = TRUE)
    }else
    {
      if(input$common_name == 1){
        curr_list <- input$print_disp_names
        remap <- nodes %>% filter(feature %in% curr_list) %>% pull(`Common Name`)
        updateSelectizeInput(session, 'print_disp_names', choices = nodes %>% pull(`Common Name`),  selected = remap, server = TRUE)
      }else{
        curr_list <- input$print_disp_names
        remap <- nodes %>% filter(`Common Name` %in% curr_list) %>% pull(feature)
        updateSelectizeInput(session, 'print_disp_names', choices = nodes %>% pull(feature), selected = remap, server = TRUE)
      }
    }
  })
  
  
  observeEvent(input$search_additional, {
    if(length(input$search_additional) == 0 ){
      if(length(input$gl) == 0){
        gene_list(NULL)
      }else{
        temp <- input$gl
        g <- sapply(temp, function(x) 
          if(x %in% genename_map$common_name){
            x <- genename_map$feature_name[which(x == genename_map$common_name)]
          }else{
            x <- x
          })
        gene_list(g)
      }
    }
    if(length(gene_list()> 0)){
      subNet <- geneListSubgraph(Net, Module, gene_list(), input$search_additional)
      Nodes <- subNet %N>% as_tibble()
      Nodes_gene_list <- Nodes %>% filter(feature %in% gene_list()) 
      if(input$common_name == 1){
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`),  selected = Nodes_gene_list %>% pull(`Common Name`), server = TRUE)
      }else{
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(feature), selected = Nodes_gene_list %>% pull(feature), server = TRUE)
      }
      sub_net(subNet)
    }
  }, ignoreNULL = FALSE)
  
  
#observeEvent(input$search_additional, {
#    if("stein" %in% input$search_addition){	
#	id <- showNotification("Generating Stiener Tree...", duration = NULL, closeButton = FALSE)
#	on.exit(removeNotification(id), add = TRUE)
#	st <- buildSteinerTrees(Net, gene_list())
#	st <- st %E>%
#	mutate(is_steiner = TRUE) %N>%
#	mutate(is_steiner = TRUE)
#	sub_net(graph_join(sub_net(), st) %>%
#             mutate(color_code = if_else(is_steiner, "#F9B6AF", "#BEBEBE")))
#	gene_list(sub_net() %N>% pull(feature))
#   }
# })
  
  observeEvent(input$method, {
    sub_net(tbl_graph())
  })
  
  observeEvent(sub_net(), {
    disp_nodes(sub_net() %N>% as_tibble())
  } )
  
  ############## Main Render ######################
  output$table  = DT :: renderDataTable({
    tabPanel("nodes", )
  })
  
  output$nodes_table <- DT::renderDataTable({
    file_name <- paste('node_table', ifelse(str_length(input$file_name) > 0, input$file_name, 'file') , sep ="_")
    S <-sub_net()
    if(isempty(S %N>% as_tibble())){
    }else{
      S_tables <- graph2NodeEdgeTables(S)
        S_nodes <- prepNodeTable(S_tables[[1]], 1) #CC: Set the # of GO-terms to 1, and removed the ability pick from the UI. To make it simpler
      DT::datatable(S_nodes, escape = FALSE, 
                    extensions = 'Buttons', options = list(
                      dom = 'Blfrtip',
                      title = paste('node_table', input$file_name, sep ="_"),
                      buttons = 
                        list('copy', 'print', list(
                          extend = 'collection',
                          buttons = list(list(extend = 'csv', filename = file_name),
                                         list(extend = 'excel', filename = file_name),
                                         list(extend = 'pdf', filename = file_name)),
                          text = 'Download')),
                      lengthMenu = list(c(10,50, 100, -1), 
                                        c('10', '30', '50', 'All')),
                      paging = T))
    }
  })
  
  output$module_table <- DT::renderDataTable({
    file_name <- paste('module_table', ifelse(str_length(input$file_name) > 0, input$file_name, 'file') , sep ="_")
    S <-sub_net()
    if(isempty(S %N>% as_tibble())){
    }else{
      S_tables <- graph2NodeEdgeTables(S)
      curr_modules <- unique(c(S_tables[[1]] %>% pull(module), unlist(S_tables[[1]] %>% pull(enriched_modules))))
      Module <- computeEnrichment(Module, S_tables[[1]] %>% pull(feature), 
                                  length(Net %N>% pull(feature)))
      ModT <- prepModuleTable(Module %>% filter(module %in% curr_modules), input$method)
      DT::datatable(ModT, escape = FALSE,
                    extensions = 'Buttons', options = list(
                      dom = 'Blfrtip',
                      title = paste('module_table', input$file_name, sep ="_"),
                      buttons = 
                        list('copy', 'print', list(
                          extend = 'collection',
                          buttons = list(list(extend = 'csv', filename = file_name),
                                         list(extend = 'excel', filename = file_name),
                                         list(extend = 'pdf', filename = file_name)),
                          text = 'Download')),
                      lengthMenu = list(c(10,50, 100, -1), 
                                        c('10', '30', '50', 'All')),
                      paging = T))
    }
  })
  
  # output$network <- renderForceNetwork({
  #   S<-sub_net()
  #   if(isempty(S %N>% as_tibble())){
  #   }else{
  #     S_tables <- graph2NodeEdgeTables(S)
  #     S_nodes <- S_tables[[1]]
  #     S_edges <- S_tables[[2]]
  #     S_edges <- S_edges %>% add_row(from = 0, to = 0, weight = 0)
  #     
  #     if(input$common_name == 1){
  #       names <- "Common Name"
  #     }else if(input$common_name == 2){
  #       names <- "feature"
  #     }
  #     
  #     if(is.null(input$show_gene_names)){
  #       op <- 0
  #       fs <- 40 
  #     }else{
  #       op <- .75
  #       fs <- 25
  #     }
  #     
  #     if(input$group_by == "regulator"){
  #       colorScale = JS('color=d3.scaleOrdinal([`#fb8072`, `#80b1d3`]), color.domain(["src","tar"])');
  #     }else if(input$group_by == "module"){
  #       num_mods <- length(setdiff(unique(S_nodes %>% pull(module)), -9999))
  #       max_colors <- palettes$maxcolors[which(palettes$pal == "Pastel1")]
  #       pal <- brewer.pal(max(3, min(max_colors, num_mods)), "Pastel1")
  #       if(max_colors < num_mods){
  #         pal <- extend_palette <- colorRampPalette(pal)(num_mods)
  #       }
  #       colorScale = JS(paste0('color=d3.scaleOrdinal([ `#BEBEBE`, ', paste(sprintf('`%s`', pal), collapse = ', '), ']), color.domain([-9999])'))
  #     }else if(input$group_by == "geneSuper"){
  #       num_supers <- length(setdiff(unique(S_nodes %>% pull(geneSuper)), "Unlabeled"))
  #       max_colors <- palettes$maxcolors[which(palettes$pal == "Pastel1")]
  #       pal <- brewer.pal(max(3, min(max_colors, num_supers)), "Pastel1")
  #       if(max_colors < num_supers){
  #         pal <- extend_palette <- colorRampPalette(pal)(num_supers)
  #       }
  #       colorScale = JS(paste0('color=d3.scaleOrdinal([`#BEBEBE`, ', paste(sprintf('`%s`', pal), collapse = ', '),"]), color.domain(['Unlabeled'])"))
  #     }
  #     
  #     if(input$method == "diff"){
  #       S_nodes <- S_nodes %>% mutate(size = rescale(score, to = c(4, 16)))
  #       forceNetwork(Links = S_edges, Nodes = S_nodes,
  #                    Source = "from", Target = "to",
  #                    Value = "weight", NodeID = names,
  #                    Group = input$group_by, Nodesize = "size", opacity = 1, opacityNoHover = op, colourScale = colorScale, 
  #                    zoom = TRUE, fontSize=fs, radiusCalculation = JS("d.nodesize"),
  #                    charge = -10) 
  #     }
  #     else{forceNetwork(Links = S_edges, Nodes = S_nodes,
  #                       Source = "from", Target = "to",
  #                       Value = "weight", NodeID = names,
  #                       Group = input$group_by, linkColour=S_edges$color_code,
  #                       opacity = 1, opacityNoHover = op, zoom = TRUE, fontSize=fs, colourScale = colorScale, 
  #                       charge = -10)
  #     }
  #   }
  # })
  
  #######  Printer Setup 
  #observeEvent(input$openPrinter, {
  #	Nodes <- disp_nodes()
  #	if(nrow(Nodes) == 0 ){
  #		showNotification("Nothing to display.")
  #		toggleModal(session, modalId ="modalPrinter", toggle = 'close')
  #	}else{
  #		if(input$common_name == 1){
  #			updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`), server = TRUE)
  #			}else{
  #			updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(feature), server = TRUE)
  #		}
  #		num_mods <- length(Nodes %>% pull(module))
  #updateSelectizeInput(session, 'print_node_pal', choices = palettes %>% pull(pal), server = TRUE)
  #	} 
  #})
  
  
  output$print_net <- renderPlot({
    subNet <- sub_net() 
    if(is_empty(subNet)){
      gg<- ggplot() + 
        theme(panel.background = element_rect(fill="white", colour = "white")) +
        geom_text(label = "no subgraph selected.")
    }else{
      
      subNet <- subNet %N>% mutate(component = group_components()) 
      keep_component <- subNet %N>% as_tibble() %>% 
        group_by(component) %>% 
        summarise(count  = n()) %>% 
        filter(count >= input$print_min_genes) %>% 
        pull(component)
      subNet <- subNet %N>% filter(component %in% keep_component )
      
      
      if(!is.null(input$print_disp_names)){
        if(input$common_name == 1){
          subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% input$print_disp_names, `Common Name`, NA))
        }else{
          subNet <-subNet %N>% mutate(display_name = ifelse(`feature` %in% input$print_disp_names, `feature`, NA))
        }
      }else{
        subNet <-subNet %N>% mutate(display_name = NA)
      }
      subNet <- subNet %N>% mutate(module = as.character(module)) %>% mutate(module = str_replace(module, '-9999', 'Unlabeled'))
      
      #edge_color_by <- ifelse('stein' %in% input$search_additional, 'is_steiner', NA)
      node_size_by <- ifelse(input$method =='diff', 'score', NA) 
      
      gg_out_plot(
        makeSubNetGraph(subNet, names_in_nodes = input$print_name_bool,
                        edge_color_by = input$edge_color_by, edge_color_palette = input$edge_color_palette, 
                        node_color_by = input$print_group_by, node_color_palette = input$print_node_pal, 
                        node_size_by = node_size_by, max_node_size = input$print_max_node_size, 
                        layout = input$print_layout, focus_nodes = list(), 
                        font_size = input$print_font_size, 
                        nudge_y = input$print_nudge_y, text_angle = input$print_text_angle, show_legend = TRUE,
                        expand_x = input$print_expand_x, expand_y = input$print_expand_y, color_scale_limits = input$edge_color_range_reg, legend_font_size = input$legend_font_size)
      )
      gg_out_plot()
    }
  })
  
  
  output$expression_heatmap <- renderPlotly({
    subNet <- sub_net() 
    if(is_empty(subNet)){
      gg<- ggplot() + 
        theme(panel.background = element_rect(fill="white", colour = "white")) +
        geom_text(label = "no subgraph selected.")
    }else{
      gg_out_heatmap(
        makeSubgraphHeatmap(subNet,
          display_name = input$common_name_heatmap,
          font_size = input$print_font_size,
          tfa_color_palette = input$tfa_palette_heatmap, 
          expression_color_palette = input$expression_palette_heatmap,
          scale_edge_color  = input$edge_color_range_heatmap, 
          scale_expression_colors = input$expression_range_heatmap, 
          scale_tfa_colors = input$tfa_range_heatmap, 
          figure_font_size = input$Font_size_heatmap 
        )
      )
      makeSubgraphHeatmapDynamic(subNet,
                          display_name = input$common_name_heatmap,
                          font_size = input$print_font_size,
                          tfa_color_palette = input$tfa_palette_heatmap, 
                          expression_color_palette = input$expression_palette_heatmap,
                          scale_edge_color  = input$edge_color_range_heatmap, 
                          scale_expression_colors = input$expression_range_heatmap, 
                          scale_tfa_colors = input$tfa_range_heatmap, 
                          figure_font_size = input$Font_size_heatmap 
      )
    }
  })
  
  
  
  
  
  output$saveFig <- downloadHandler(
    filename = function(){ifelse(str_length(input$print_file_name) > 0, 
                                 paste(input$print_file_name, '.', input$print_file_type, sep = ''), 
                                 paste('file', '.', input$print_file_type, sep = ''))
    },
    content = function(file){
      ggsave(file,gg_out_plot(), width = input$print_image_width, height = input$print_image_height,units = 'in')
    })
  
  
  output$saveFig2 <- downloadHandler(
    filename = function(){ifelse(str_length(input$print_file_name) > 0, 
                                 paste(input$print_file_name, '.', input$print_file_type, sep = ''), 
                                 paste('file', '.', input$print_file_type, sep = ''))
    },
    content = function(file){
      ggsave(file, gg_out_heatmap(), width = input$print_image_width_heatmap, height = input$print_image_height_heatmap,units = 'in')
    })
  
  
  
  observeEvent(input$plot_click, {
    Nodes <- disp_nodes()
    gg_out<-gg_out_plot()
    gg_data<-tibble(gg_out$data)
    subNet <- sub_net()
    if(!is_empty(subNet)){
      if(input$common_name == 1){
        gg_name  <- nearPoints(gg_data, input$plot_click, threshold = 35, maxpoints = 1) %>% pull(`Common Name`)
      }else{
        gg_name <- nearPoints(gg_data, input$plot_click, threshold = 35, maxpoints = 1) %>% pull(feature)
      }
      
      if(!is_empty(gg_name)){
        if(gg_name %in% input$print_disp_names){
          updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`), server = TRUE, selected = setdiff(input$print_disp_names, gg_name))
        }else{
          if(input$common_name == 1){
            updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`), server = TRUE, selected = c(unlist(input$print_disp_names), gg_name))
          }else{
            updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(feature), server = TRUE, selected = c(unlist(input$print_disp_names), gg_name))
          }
        }
      }
    }
  })
  
  
  ############### Save Features #####################################
  output$save_node_info <- downloadHandler(
    filename = function() {
      paste(node_name_info(), '.txt', sep='')
    },
    content = function(con) {
      write(str_replace_all(printNodeInfo(Net, node_name_info()), '<br/>', '\n'), con)
    }
  )
  
  output$save_module_info <- downloadHandler(
    filename = function() {
      paste("module_",module_id_info(), '.txt', sep='')
    },
    content = function(con) {
      text_info <- printModuleInfo(Module, module_id_info(), list())
      text_info <- str_replace_all(text_info, '<br/>', '\n')
      text_info <- str_replace_all(text_info, '&emsp', '\t')
      write(text_info, con)
    }
  )
  
  output$save_all_module_info <- downloadHandler(
    filename = function() {
      if(input$method =="go_term"){
        paste(input$go_term, "enriched_modules.txt", sep='')  
      }
      else if(input$method == "gene"){
        paste(input$gene, "_modules.txt", sep='')
      }
      else if(input$method == "list"){
        file <- input$cell_list_file
        fp <- file$name
        paste(file_path_sans_ext(basename(fp)),"_modules.txt", sep='')  
      }
    },
    content = function(con) {
      if(input$method == "list"){
        text_info <- printAllModuleInfo(sub_net(), Module, gene_list(), genes)
      }else{
        text_info <- printAllModuleInfo(sub_net(), Module, list(), genes)
      }
      text_info <- str_replace_all(text_info, '<br/>', '\n')
      text_info <- str_replace_all(text_info, '&emsp;', '\t')
      write(text_info, con)
    }
  )
  
  ### Initialize 
  disp_names <- c('cyp51A', 'erG25B', 'hyd1', 'srbA', 'srbB', 'erG3', 'erG25', 'fhpA', 'erG1', 
                       'hem13', 'niiA', 'AFUA_5G06120_nca', 'AFUA_3G12190', 'bna4', 'srb5', 'hem14', 'exG4', 
                       'erG3A', 'pre4', 'AFUA_7G04740', 'AFUA_6G02180')
  init <<- TRUE
  ############### observe Events ##########################################
  #observeEvent(input$node_name, {
  #  if(any(input$node_name == genename_map$common_name)){
  #    idx <- which(input$node_name == gene_map)
  #    name <- genename_map$feature_name[idx]
  #  }else{
  #    name <- input$node_name
  #  }
  #  node_name_info(name)
  #  if(input$method!="module"){
  #    module_id_info(getModuleID(Net, node_name_info()))
  #  }
  #})
  
  #output$node_info <- renderUI({
  #    HTML(printNodeInfo(Net, node_name_info()))
  #})
  
  
  #output$module_info <-renderUI({
  #  if(input$method == "list"){
  #    gl <- gene_list()
  #  }else{
  #    gl <- list()
  #  }
  #  text <- printModuleInfo(Module, module_id_info(), gl, genes)
  #  HTML(text)
  #})

########################################################################


showModal(
  modalDialog(
    title = "Welcome to GRAsp",
    "If this is your first time using GRAsp, please read the documentation in the About tab.",
    footer = tagList(
      actionButton("go_to_about", "Go to About", class = "btn-primary"),
      modalButton("Close")  # This adds a close button
    ),
    easyClose = TRUE
  )
)

# Observe button click and switch to About tab
observeEvent(input$go_to_about, {
  updateNavbarPage(session = getDefaultReactiveDomain(), "navbar", selected = "About")  # Switch tab
  removeModal()  # Close the popup
})

}
shinyApp(ui, server)


