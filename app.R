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


## Load in aux functions and data ----
source('aux_functions.R') 
source('printerFunction.R')
source('heatmap_printer.R')
source('heatmap_printer_dynamic.R')
source('ui_items.R')

## Initialize global variables ----
all_gene_names <- unique(c(genes, genename_map$common_name))
palettes_nodes<- tibble(rownames_to_column(brewer.pal.info, var = 'pal')) 
palettes_edges <- tibble(rownames_to_column(brewer.pal.info, var = 'pal')) %>% 
  filter(category %in% c("div", "seq"))
igraph_layout <- c('Fruchterman-Reingold'='nicely', 'Davidson-Harel'='dh', 'Kamada-Kawai'='kk', 'Large graph layout'= 'lgl') #'Force directed' = 'drl')

mean_names <- names((Net %N>% pull(expression))[[1]]) 
cluster <- unique(str_split_i(mean_names, '-', 1))
type <- unique(str_split_i(mean_names, '-', 2))

# ui  -----
ui <- navbarPage(title,
                 id = 'navbar',
                 theme = shinytheme("flatly"), #CC: Set this as the main theme
                 
                 tabPanel("Visualize",
                          ## Left bar setting ----- 
                          fluidRow(
                            column(2,
                                   ### Method selection ------
                                   method_Picker(),
                                   conditionalPanel(
                                     condition = "input.method == 'module'",
                                     module_Picker(module_ids = module_ids)
                                   ),
                                   ### GO Selection -------
                                   conditionalPanel(
                                     condition = "input.method == 'go_term'",
                                     pickerInput(inputId = "go_term", label = "GO terms", choices = c("", unlist(sort(enriched_go_terms))))
                                   ),
                                   
                                   ### Gene List selection ------
                                   conditionalPanel(
                                     condition = "input.method == 'list'",
                                     geneList_Picker_Tag(),
                                     geneList_Picker(),
                                     geneList_File_Tag(), 
                                     geneList_File(),
                                     geneList_checkBox()
                                   ),
                                   
                                   ### Diffusion Selection ------
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
                            ## Visualization Block ########
                            column(10,
                                   tabsetPanel(id = "displayType", type = "tabs",
                                               tabPanel("Network Plot", plotOutput("print_net", click = 'plot_click', height = '1000px')),
                                               tabPanel("Expression Heatmaps", plotlyOutput("expression_heatmap", height = "1000px")),
                                               tabPanel("Gene Table", DT::dataTableOutput("nodes_table")),
                                               tabPanel("Module Table", DT::dataTableOutput("module_table"))
                                               
                                   )
                            )
                          ),
                          fluidRow(
                            ## Network plot support options  -----
                            conditionalPanel(
                              condition ="input.displayType == 'Network Plot'",
                              column(2),
                              ### Graph visualization setting -------
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
                              
                              ### Sample selection 
                              column(2,
                                     networkViz_sampleSelect_Tag(),
                                     networkViz_nodeExpGlobal_Checkbox(),
                                     conditionalPanel(condition = "input.global == false", 
                                       networkViz_dispSample_Select(cluster)
                                     )
                              ),
                              ### Node visualization settings --------
                              column(2,
                                     networkViz_nodeColor_Radio(), 
                                     conditionalPanel(condition = "input.print_group_by == 'exp'",
                                        networkViz_expColor_Slider(limits = default_expression_range)
                                     ),
                                     networkViz_nodeColor_Select(palettes_nodes),
                                     networkViz_nodeSize_Slider(),
                                     networkViz_nodeFontSize_Slider()
                              ),
                              
                              ### Edge visualzation setting ------
                              column(2,
                                     # Removing this and having correlation scale width
                                     # networkViz_edgeColor_Radio(),
                                     # conditionalPanel(
                                     #   condition = "input.edge_color_by == 'Reg_weight'",
                                     #   networkViz_edgeRangeReg_Slider()
                                     # ),
                                     # conditionalPanel(
                                     #   condition = "input.edge_color_by == 'Correlation'",
                                     #   networkViz_edgeColorRangeCorr_Slider()
                                     # ),
                                     networkViz_edgeRangeReg_Slider(),
                                     networkViz_edgeRangeAbsCorr_Slider(),
                                     #networkViz_edgeWidth_Slider(), 
                                     networkViz_edgePalette_Select(palettes_edges)
                               ),
                              
                              ### Save Setting -------
                              column(2,
                                     networkViz_expandX_Slider(),
                                     networkViz_expandY_Slider(),
                                     networkViz_legendSize_Slider(),
                                     networkViz_imageSaveWidth_Slider(),
                                     networkViz_imageSaveHeight_Slider(),
                                     networkViz_imageName_Text(),
                                     networkViz_imageType_Select(),
                                     networkViz_Download(), 
                              )
                            ),
                            ## Heatmap support option tab  -------
                            conditionalPanel(
                              condition ="input.displayType == 'Expression Heatmaps'",
                              column(2),
                              column(2,
                                     heatmapViz_nameFormat_Radio()
                              ),
                              column(2,
                                     #heatmapViz_TFAPalette_Select(),
                                     #heatmapViz_TFARange_Slider(), 
                                     heatmapViz_expPalette_Select(palettes_edges), 
                                     heatmapViz_expRange_Slider(),
                                     heatmapViz_direction_Radio()
                              ),
                              
                              ## commenting out heatmap options relating to TFA
                              ## and network visualization in heatmap panel. 
                              # column(2,
                              #        heatmapGraphViz_edgeColor_Radio(),
                              #        conditionalPanel(
                              #          condition = "input.edge_color_by_heatmap == 'Reg_weight'",
                              #          heatmapGraphViz_edgeRangeReg_Slider()
                              #        ),
                              #        conditionalPanel(
                              #          condition = "input.edge_color_by_heatmap == 'Correlation'",
                              #          heatmapGraphViz_edgeRangeCorr_Slider()
                              #        ),
                              #        heatmapGraphViz_edgePalette_Select(palettes_edges)
                              # ),
                              column(2,
                                     heatmapViz_fontSize_Slider(),
                                     heatmapViz_imageWidth_Slider(),
                                     heatmapViz_imageHeight_Slider(), 
                              ),
                              column(2, 
                                     heatmapViz_fileName_Text(),
                                     heatmapViz_fileType_Select(),
                                     heatmapViz_Download()
                              )
                            )
                          )
                 ),
                 
                 ## Additional tab support -----
                 tabPanel("About",
                          fluidPage(htmltools::tags$iframe(src = "Help.html", width = '100%', height = 1000, style = "border:none;"))
                 ),
                 tabPanel("Contact", #CC: connected a googleforms that allows users to ask questions or give suggestions. Not sure if there was a better way.
                          tags$iframe(
                            src="https://docs.google.com/forms/d/e/1FAIpQLSdmkuOZk7DCRHWBw5cOG028fvMjNO9yZL3L3FZrDhtyppN2pQ/viewform?embedded=true",
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
# Server Functions ########################################
server <- function(input, output, session) {
  ## Initialize -----
  disp_names <- c('cyp51A', 'erG25B', 'hyd1', 'srbA', 'srbB', 'erG3', 'erG25', 'fhpA', 'erG1', 
                  'hem13', 'niiA', 'AFUA_5G06120_nca', 'AFUA_3G12190', 'bna4', 'srb5', 'hem14', 'exG4', 
                  'erG3A', 'pre4', 'AFUA_7G04740', 'AFUA_6G02180')
  init <- TRUE
  
  
  ## Popup ----
  showModal(
    modalDialog(
      title = "Welcome to MERLIN-VIZ",
      "If this is your first time using MERLIN-VIZ, please read the documentation in the About tab.",
      footer = tagList(
        actionButton("go_to_about", "Go to About", class = "btn-primary"),
        modalButton("Close")  # This adds a close button
      ),
      easyClose = TRUE
    )
  )
  
  ### Observe button click and switch to About tab ----
  observeEvent(input$go_to_about, {
    updateNavbarPage(session = getDefaultReactiveDomain(), "navbar", selected = "About")  # Switch tab
    removeModal()  # Close the popup
  })
  
  
  ## Dynamic selectize update ----
  updateSelectizeInput(session, 'gene', choices = all_gene_names, server = TRUE)
  updateSelectizeInput(session, 'gl', choices = all_gene_names, selected = default_gene, server = TRUE)
  
  
  ## Dynamic Variables ------
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
  
  
  ## File IO ----
  ### Gene list file path ----
  file_path <- reactive({
    file <- input$cell_list_file
    file$datapath 
  })
  
  ### diffusion file path ----
  diff_file_path <- reactive({
    file<-input$diff_list_file
    req(file)
    file$datapath
  })
  
  ### Reset diffusion button ----
  observeEvent(input$refresh_diff, {
    percentile(input$percentile)
    lambda(input$kernel)
  })
  
  ### Diffusion min neighbors button ---- 
  observeEvent(input$min_neigh, {
    min_neigh(input$min_neigh)
    render_diff(TRUE)
  })
  
  ### Diffusion regulator number ---- 
  observeEvent(input$disp_regs, {
    disp_regs(input$disp_regs)
    render_diff(TRUE)
  })
  
  ### Re-render diffusion plots ----
  observeEvent(render_diff(), {
    nodes <- diff_nodes()
    if(!is.null(nodes$score)){
      Net <- left_join(Net, nodes)
      sub_net(diffScoreSubgraph(Net, min_neigh(), disp_regs()))
    }
    render_diff(FALSE)
    #print(render_diff())
  })
  
  
  ## Update steps for all subnets -----
  ### Gene list update -----
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
        }else{
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
  
  ### Go term update ----
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
   
  ### Module search update ----
  observeEvent(input$module_id,{
    if(input$module_id  == ""){
      sub_net(tbl_graph())
    }else{
      module <- input$module_id
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
  
  ### Diffusion file update -----
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
  
  ### Gene list from file update  -----
  observeEvent(file_path(), {
    fp <- file_path()
    sub_net(gene_list(read_csv(file = fp, col_names=FALSE) %>% pull(X1)))
  })
  
  ### Gene list subgraph generation -----
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
  
  ### Search additional update ----- 
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

  ### Clear network update ----  
  observeEvent(input$method, {
    sub_net(tbl_graph())
  })
  
  ### Clear network display gene name  update ----
  observeEvent(sub_net(), {
    disp_nodes(sub_net() %N>% as_tibble())
  } )
  
  ## Main Renders -----
  ### Node table render ----
  output$nodes_table <- DT::renderDataTable({
    file_name <- paste('node_table', ifelse(str_length(input$file_name) > 0, input$file_name, 'file') , sep ="_")
    S <-sub_net()
    if(isempty(S %N>% as_tibble())){
    }else{
      S_tables <- graph2NodeEdgeTables(S)
        S_nodes <- prepNodeTable(S_tables[[1]], 1) 
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
  
  ### Module table render ---- 
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

  
  ### Render network visualization ----
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
      
      ## Set up node feature selection block
      if(input$print_group_by== "exp"){
        if(input$global){
          node_color_by <- "mean_expression"
        }else
        {
          cluster_sample_string <- input$disp_sample #paste(input$disp_cluster, input$disp_type, sep = '-')
          subNet <- subNet %N>% mutate(exp = map_dbl(expression, ~ .x[cluster_sample_string]))
          node_color_by <- "exp"
        }
      }else{
        node_color_by <- input$print_group_by
      }
      
      ## Set up block for feature selection for edges. 
      if(input$global){
        edge_color_by <- 'Reg_weight'
        edge_width_by <- 'Correlation'
      }else{
        cluster_sample_string <- input$disp_sample 
        edge_color_by <- paste(cluster_sample_string, 'Reg_weight', sep = '_')
        edge_width_by <- paste(cluster_sample_string, 'Correlation', sep = '_')
      }
      
      if(!is.null(input$print_disp_names)){
        if(input$common_name == 1){
          subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% input$print_disp_names, `Common Name`, NA))
        }else{
          subNet <-subNet %N>% mutate(display_name = ifelse(`feature` %in% input$print_disp_names, `feature`, NA))
        }
      }else{
        subNet <-subNet %N>% mutate(display_name = NA)
      }
      subNet <- subNet %N>% 
        mutate(module = as.character(module)) %>% 
        mutate(module = str_replace(module, '-9999', 'Unlabeled')) %>% 
        mutate(subnet_degree = centrality_degree(mode = "all", weights = NULL))
      
      node_size_by <- ifelse(input$method =='diff', 'score', 'subnet_degree') 
      
      # if(input$edge_color_by == "Reg_weight"){
      #   edge_color_range = input$edge_color_range_reg
      # }else{
      #   edge_color_range = input$edge_color_range_corr
      # }
      edge_color_range = input$edge_color_range_reg
      
      gg_out_plot(
        makeSubNetGraph(subNet, names_in_nodes = input$print_name_bool, node_color_by = node_color_by, 
                        edge_color_by = edge_color_by, edge_color_palette = input$edge_color_palette, 
                        edge_width_by = edge_width_by, max_edge_width = input$edge_range_abs_corr, 
                        node_color_palette = input$print_node_pal, 
                        node_size_by = node_size_by, max_node_size = input$print_max_node_size, 
                        layout = input$print_layout, focus_nodes = list(), 
                        font_size = input$print_font_size, 
                        nudge_y = input$print_nudge_y, text_angle = input$print_text_angle, show_legend = TRUE,
                        expand_x = input$print_expand_x, expand_y = input$print_expand_y, color_scale_limits = edge_color_range, 
                        node_scale_limits = input$exp_color_scale,
                        legend_font_size = input$legend_font_size, direction = 1)
      )
      gg_out_plot()
    }
  })
  
  ### Render expression heatmap ---- 
  output$expression_heatmap <- renderPlotly({
    subNet <- sub_net() 
    d <- input$heatmap_direction 
    if(is_empty(subNet)){
      gg<- ggplot() + 
        theme(panel.background = element_rect(fill="white", colour = "white")) +
        geom_text(label = "no subgraph selected.")
    }else{
      print(input$tfa_color_range_heatmap)
      gg_out_heatmap(
        #### Static for output figures ----
        makeSubgraphHeatmap(subNet,
                            display_name = input$common_name_heatmap,
                            #edge_color_by = input$edge_color_by_heatmap, 
                            #edge_color_palette = input$edge_color_palette_heatmap, 
                            font_size = input$print_font_size,
                            direction = input$heatmap_direction, 
                            tfa_color_palette = input$tfa_palette_heatmap, 
                            expression_color_palette = input$expression_palette_heatmap,
                            scale_edge_color  = input$edge_color_range_heatmap, 
                            scale_expression_colors = input$expression_range_heatmap, 
                            scale_tfa_colors = input$tfa_range_heatmap, 
                            figure_font_size = input$Font_size_heatmap 
        )
      )
      
      #### Dynamic within app ----
      makeSubgraphHeatmapDynamic(subNet,
                                 display_name = input$common_name_heatmap,
                                 #edge_color_by = input$edge_color_by_heatmap, 
                                 #edge_color_palette = input$edge_color_palette_heatmap, 
                                 font_size = input$print_font_size,
                                 direction = as.integer(input$heatmap_direction),
                                 tfa_color_palette = input$tfa_palette_heatmap, 
                                 expression_color_palette = input$expression_palette_heatmap,
                                 scale_edge_color  = input$edge_color_range_heatmap, 
                                 scale_expression_colors = input$expression_range_heatmap, 
                                 scale_tfa_colors = input$tfa_range_heatmap, 
                                 figure_font_size = input$Font_size_heatmap 
      )
    }
  })
  
  
  ## Download handlers -----
  ### Network download handler -----
  output$saveFig <- downloadHandler(
    filename = function() {
      ext <- input$print_file_type
      name <- if (str_length(input$print_file_name) > 0) input$print_file_name else "file"
      paste0(name, ".", ext)
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = gg_out_plot(),
        width = input$print_image_width,
        height = input$print_image_height,
        units = "in",
        device = input$print_file_type
      )
    }
  )
  
  
  ### heatmap download handler -----
  output$saveFig2 <- downloadHandler(
    filename = function() {
      ext <- input$print_file_type
      name <- if (str_length(input$print_file_name) > 0) input$print_file_name else "file"
      paste0(name, ".", ext)
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = gg_out_heatmap(),
        width = input$print_image_width_heatmap,
        height = input$print_image_height_heatmap,
        units = "in",
        device = input$print_file_type2
      )
    }
  )
  
  
  
  ##  Observe event for name update ----
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
}
shinyApp(ui, server)


