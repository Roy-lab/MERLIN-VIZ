library(shiny)
library(tools)
library(tidyverse)
library(networkD3)
library(Matrix)
library(scales)
library(DT)
library(webshot)
library(htmlwidgets)
library(RColorBrewer)
library(shinyBS)

## Load in aux functions and data
source('aux_functions.R')
source('printerFunction.R')
#MyClickScript <- 'Shiny.setInputValue("save_module", module_name)'

all_gene_names <- unique(c(genes, genename_map$common_name))
palettes <- tibble(rownames_to_column(brewer.pal.info, var = 'pal'))
igraph_layout <- c('Fruchterman-Reingold'='nicely', 'Davidson-Harel'='dh', 'Kamada-Kawai'='kk', 'Large graph layout'= 'lgl', 'Force directed' = 'drl')

################################# ui ###########################################
ui <- navbarPage("GRAsp",
  tabPanel("Visualize",
    fluidRow(
      column(2, 
             radioButtons("common_name", h4("Gene Name Format"),
                          choices = list("Common" = 1, "Standard" = 2)),
             radioButtons("group_by", h4("Node Color By"), 
                          choices = list ("Module" = 'module', "Regulator" = "regulator", "Gene Family"= 'geneSuper')),
             selectInput(inputId = "method", h4("Search Method"),
                         c("", Modules = "module", 'GO Term' = "go_term", 
                           'Gene List' = 'list', "Node Diffusion" = "diff")
             ),
             conditionalPanel(
               condition = "input.method == 'module'",
               selectInput(
                 inputId = "module_id", "Module ID",
                 sort(module_ids)
               )
             ),
             conditionalPanel(
               condition = "input.method == 'go_term'",
               selectInput(inputId = "go_term", label = "GO terms", choices = sort(enriched_go_terms))
             ),
             conditionalPanel(
               condition = "input.method == 'list'",
               radioButtons("list_type", h4("Select List Type"), choices = list("Select" = 1, "File" = 2)),
               conditionalPanel(
                 condition = "input.list_type == 1", 
                 selectInput(inputId = "gl", label = "Genes", choices = NULL, multiple = TRUE, selectize = TRUE)
               ),
               conditionalPanel(
                 condition = "input.list_type == 2", 
                 fileInput(inputId = "cell_list_file", "Cell List File")
               ),
               checkboxGroupInput(inputId = "search_additional", 
                                  label = h4("Include Additional Genes"),
                                  c("Steiner Tree" = "stein", "Modules" = "mod", "Neighbors" = "neigh")
               )
             ),
             conditionalPanel(
               condition ="input.method == 'diff'",
               fileInput(inputId = "diff_list_file", "Diffusion Score File"),
               conditionalPanel(
                 condition = "output.diffFileUploaded",
                 selectInput(inputId = "kernel", "Lambda", c(1, 10, 100, 1000)),
                 numericInput(inputId = "min_neigh", label = "Minimum Number of Targets", value = 5),
                 numericInput(inputId = "disp_regs", label = "Number of Regulators to Display", value = 5),
                 #sliderInput(inputId = "percentile", "Node Percentile", min = 90, max = 100,
                 #          value = 95, step = .1),
                 actionButton("refresh_diff", "Refresh")
               )
             ),
             checkboxGroupInput(inputId = "show_gene_names", label = h4("Show Gene Names"), 
                                c("Show Names" = "sh")),
             numericInput(inputId = "disp_go", "Number of GO in Table", value = 5),
             
             actionButton("openPrinter", "Save Network")
             
      ),
      column(10,forceNetworkOutput("network",height = "750px")),
      bsModal("modalPrinter", "Print Image", "openPrinter", size = "large", 
              fluidRow(
                plotOutput(outputId = "print_net", click = 'plot_click'),
                verbatimTextOutput("node_name")
              ),
              fluidRow(
                column(4,
                      selectInput(inputId = 'print_layout', label = 'Node layout', choices = igraph_layout, multiple = FALSE),
                      numericInput(inputId = 'print_min_genes', label = "Minimum number of genes in component", value = 0, min = 0),
                      selectInput(inputId = 'print_disp_names', label = 'Names to display', choices = NULL, multiple = TRUE, selectize = TRUE),
                      checkboxInput(inputId = 'print_name_bool', label = "Name in nodes", value = TRUE),
                      conditionalPanel(
                        condition = "input.print_name_bool == 0", 
                        numericInput(inputId = 'print_nudge_y', label = "Nudge labels", value = 0.3, min = 0, step = 0.1),
                        sliderInput(inputId = 'print_text_angle', label = "Text angle", value = 0, min = -90, max = 90, step = 15)
                      )
                  ),
                  column(4,
                      radioButtons(inputId = "print_group_by", label = "Node color by", 
                                  choices = list ("Module" = 'module', "Regulator" = "regulator", "Gene Family"= 'geneSuper'), selected =  'module'),
                      selectInput(inputId = 'print_node_pal', label = 'Color pallette', choices = palettes$pal, multiple = FALSE, selected = 'Pastel1'),
                      numericInput(inputId = 'print_max_node_size', label = "Node size", value = 8, min = 1),
                      numericInput(inputId = 'print_font_size', label = "Font size", value = 8, min = 1)
                  ),
                  column(4,
                      numericInput(inputId = 'print_expand_x', label = "Expand X", value = 0, step = 0.1),
                      numericInput(inputId = 'print_expand_y', label = 'Expand Y', value = 0, step = 0.1),
                      numericInput(inputId = 'print_image_height', label = 'Image height (in)', value = 8, min = 1),
                      numericInput(inputId = 'print_image_width', label = 'Image width (in)', value = 8, min = 1),
                      textInput("print_file_name", "Figure name", value = "file_name"),
                      downloadButton("saveFig", "Save figure") 
                  ))
    )
    ),
    fluidRow(
      textInput("file_name", "Table File Name", value = "file_name"),
      column(12,
             tabsetPanel(
               id = 'table',
               tabPanel("nodes", DT::dataTableOutput("nodes_table")),
               tabPanel("modules", DT::dataTableOutput("module_table"))
               )
      )
    )
  ),
  tabPanel("Help",
           fluidPage(htmltools::tags$iframe(src = "Help.html", width = '100%',  height = 1000,  style = "border:none;")) 
  )
)
##################### Server Functions ########################################
server <- function(input, output, session) {
  updateSelectizeInput(session, 'gene', choices = all_gene_names, server = TRUE)
  updateSelectizeInput(session, 'gl', choices = all_gene_names, server = TRUE)
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
  disp_nodes <- reactiveVal(NULL)
  
  lambda <- reactiveVal(1)
  gene_list <- reactiveVal(NULL)
  
  
  ##### File IO ################
  file_path <- reactive({
    file <- input$cell_list_file
    print(file)
    req(file)
    file$datapath 
  })
  
  diff_file_path <- reactive({
    file<-input$diff_list_file
    req(file)
    file$datapath
  })
  
  observeEvent(input$refresh_diff, {
    percentile(input$percentile)
    print(percentile())
    lambda(input$kernel)
    print(lambda())
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
    print("IN Func")
    nodes <- diff_nodes()
    if(!is.null(nodes$score)){
      Net <- left_join(Net, nodes)
      sub_net(diffScoreSubgraph(Net, min_neigh(), disp_regs()))
    }
    render_diff(FALSE)
    print(render_diff())
  })
  
  observeEvent(input$gl,{
    print(length(Net %N>% pull(feature)))
    print(input$gl)
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
  
  observeEvent(input$go_term,{
    go_term <- input$go_term
    print(go_term)
    subgraph <- goSubgraph(Net, Module, enrich_2_module, go_term)
    sub_net(goSubgraph(Net, Module, enrich_2_module, go_term))
  })
  
  observeEvent(input$module_id,{
    module <- input$module_id
    print(module)
    module_id_info(input$module_id)
    node_name_info(NA)
    sub_net(moduleSubgraph(Net, Module, module))
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
        print("HERE:::")
        render_diff(TRUE)
  })
  
  observeEvent(file_path(), {
    fp <- file_path()
    sub_net(gene_list(read_csv(file = fp, col_names=FALSE) %>% pull(X1)))
  })
  
  observeEvent(gene_list(), {
    if(length(gene_list()> 0 )){
      sub_net(geneListSubgraph(Net, Module, gene_list(), input$search_additional))
    }
  })
  
  observeEvent(input$search_additional, {
    if(length(input$search_additional) == 0 ){
      if(length(input$gl) ==0){
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
    sub_net(geneListSubgraph(Net, Module, gene_list(), input$search_additional))
    }
  }, ignoreNULL = FALSE)
  
  
  observeEvent(input$steiner, {
    id <- showNotification("Generating Stiener Tree...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    st <- buildSteinerTrees(Net, gene_list())
    st <- st %E>%
      mutate(is_steiner = TRUE) %N>%
      mutate(is_steiner = TRUE)
    sub_net(graph_join(sub_net(), st) %>%
      mutate(color_code = if_else(is_steiner, "#F9B6AF", "#BEBEBE")))
    gene_list(sub_net() %N>% pull(feature))
  })
  
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
    S_nodes <- prepNodeTable(S_tables[[1]], input$disp_go)
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
    ModT <- prepModuleTable(Module %>% filter(module %in% curr_modules), input$method, input$disp_go)
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
  
  output$network <- renderForceNetwork({
    S<-sub_net()
    if(isempty(S %N>% as_tibble())){
    }else{
    S_tables <- graph2NodeEdgeTables(S)
    S_nodes <- S_tables[[1]]
    S_edges <- S_tables[[2]]
    S_edges <- S_edges %>% add_row(from = 0, to = 0, weight = 0)
    
    if(input$common_name == 1){
      names <- "Common Name"
    }else if(input$common_name == 2){
      names <- "feature"
    }
    
    if(is.null(input$show_gene_names)){
      op <- 0
      fs <- 40 
    }else{
      op <- .75
      fs <- 25
    }
    
    if(input$group_by == "regulator"){
      colorScale = JS('color=d3.scaleOrdinal([`#fb8072`, `#80b1d3`]), color.domain(["src","tar"])');
    }else if(input$group_by == "module"){
      num_mods <- length(setdiff(unique(S_nodes %>% pull(module)), -9999))
      max_colors <- palettes$maxcolors[which(palettes$pal == "Pastel1")]
      pal <- brewer.pal(max(3, min(max_colors, num_mods)), "Pastel1")
      if(max_colors < num_mods){
        pal <- extend_pallette <- colorRampPalette(pal)(num_mods)
      }
     colorScale = JS(paste0('color=d3.scaleOrdinal([ `#BEBEBE`, ', paste(sprintf('`%s`', pal), collapse = ', '), ']), color.domain([-9999])'))
    }else if(input$group_by == "geneSuper"){
      num_supers <- length(setdiff(unique(S_nodes %>% pull(geneSuper)), "Unlabeled"))
      max_colors <- palettes$maxcolors[which(palettes$pal == "Pastel1")]
      pal <- brewer.pal(max(3, min(max_colors, num_supers)), "Pastel1")
      if(max_colors < num_supers){
        pal <- extend_pallette <- colorRampPalette(pal)(num_supers)
      }
      colorScale = JS(paste0('color=d3.scaleOrdinal([`#BEBEBE`, ', paste(sprintf('`%s`', pal), collapse = ', '),"]), color.domain(['Unlabeled'])"))
    }

    if(input$method == "diff"){
      S_nodes <- S_nodes %>% mutate(size = rescale(score, to = c(4, 16)))
      forceNetwork(Links = S_edges, Nodes = S_nodes,
                   Source = "from", Target = "to",
                   Value = "weight", NodeID = names,
                   Group = input$group_by, Nodesize = "size", opacity = 1, opacityNoHover = op, colourScale = colorScale, 
                   zoom = TRUE, fontSize=fs, radiusCalculation = JS("d.nodesize"),
                   charge = -10) 
    }
    else{forceNetwork(Links = S_edges, Nodes = S_nodes,
                 Source = "from", Target = "to",
                 Value = "weight", NodeID = names,
                 Group = input$group_by, linkColour=S_edges$color_code,
                 opacity = 1, opacityNoHover = op, zoom = TRUE, fontSize=fs, colourScale = colorScale, 
                 charge = -10)
    }
    }
  })
  
  #######  Printer Setup 
  observeEvent(input$openPrinter, {
    Nodes <- disp_nodes()
    if(nrow(Nodes) == 0 ){
      showNotification("Nothing to display.")
      toggleModal(session, modalId ="modalPrinter", toggle = 'close')
    }else{
      if(input$common_name == 1){
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(`Common Name`), server = TRUE)
      }else{
        updateSelectizeInput(session, 'print_disp_names', choices = Nodes %>% pull(feature), server = TRUE)
      }
      num_mods <- length(Nodes %>% pull(module))
      #updateSelectizeInput(session, 'print_node_pal', choices = palettes %>% pull(pal), server = TRUE)
    } 
  })
  
  
 output$print_net <- renderPlot({
   subNet <- sub_net() 
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
   
   edge_color_by <- ifelse('stein' %in% input$search_additional, 'is_steiner', NA)
   node_size_by <- ifelse(input$method =='diff', 'score', NA) 
   
   gg_out_plot(
   makeSubNetGraph(subNet, names_in_nodes = input$print_name_bool, node_color_by = input$print_group_by, 
                                  edge_color_by = edge_color_by, node_color_palette = input$print_node_pal, 
                                  node_size_by = node_size_by, max_node_size = input$print_max_node_size, 
                                  layout = input$print_layout, focus_nodes = list(), 
                                  font_size = input$print_font_size, 
                                  nudge_y = input$print_nudge_y, text_angle = input$print_text_angle, show_legend =input$print_show_legend,
                                  expand_x = input$print_expand_x, expand_y = input$print_expand_y)
   )
   gg_out_plot()
   })
 
 
 
 
 
 output$saveFig <- downloadHandler(
   filename = function(){ifelse(str_length(input$print_file_name) > 0, paste(input$print_file_name, '.png', sep = ''), 'file.png')},
   content = function(file){
     ggsave(file,gg_out_plot(), width = input$print_image_width, height = input$print_image_height,units = 'in')
   })
 
 
 
 observeEvent(input$plot_click, {
   Nodes <- disp_nodes()
   gg_out<-gg_out_plot()
   gg_data<-tibble(gg_out$data)
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
}
########################################################################
shinyApp(ui, server)

