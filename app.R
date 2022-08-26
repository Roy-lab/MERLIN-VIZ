library(shiny)
library(tools)
library(tidyverse)
library(networkD3)
library(Matrix)
library(scales)
library(DT)
library(webshot)
library(htmlwidgets)

## Load in aux functions and data
source('asp_aux_functions.R')
#MyClickScript <- 'Shiny.setInputValue("save_module", module_name)'

all_gene_names <- unique(c(genes, genename_map$common_name))

################################# ui ###########################################
ui <- fluidPage(
    fluidRow(
      column(2, 
             radioButtons("common_name", h4("Gene Name Format"),
                          choices = list("Common" = 1, "Standard" = 2)),
             radioButtons("group_by", h4("Node Color By"), 
                          choices = list ("Module" = 'module', "Regulator" = "regulator")),
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
               selectInput(
                 inputId = "go_term", "GO terms",
                 sort(enriched_go_terms)
               )
             ),
             conditionalPanel(
               condition = "input.method == 'list'",
               radioButtons("list_type", h2("Select List Type"), choices = list("Select" = 1, "File" = 2)),
               conditionalPanel(
                 condition = "input.list_type == 1", 
                 selectInput(inputId = "gl", label = "Genes", choices = NULL, multiple = TRUE, selectize = TRUE)
               ),
               conditionalPanel(
                 condition = "input.list_type == 2", 
                 fileInput(inputId = "cell_list_file", "Cell List File")
               ),
               checkboxGroupInput(inputId = "search_additional", 
                                  label = h4("Inlude additional Genes"),
                                  c("Modules" = "mod", "Neighbors" = "neigh")
               )
             ),
             conditionalPanel(
               condition ="input.method == 'diff'",
               fileInput(inputId = "diff_list_file", "Diffusion Score File"),
               conditionalPanel(
                 condition = "output.diffFileUploaded",
                 selectInput(inputId = "kernel", "Lambda", c(1, 10, 100, 1000)),
                 sliderInput(inputId = "percentile", "Node Percentile", min = 90, max = 100,
                             value = 95, step = .1),
                 actionButton("refresh_diff", "Refresh")
               )
             ),
             checkboxGroupInput(inputId = "show_gene_names", label = h4("Show Gene Names"), 
                                c("Show Names" = "sh")),
             actionButton("steiner", "Generate Steiner Tree"),
             textInput("file_name", h4("Save File"), value = "file_name"),
             downloadButton("save_net_image", "Save Network")
             
      ),
      column(10,forceNetworkOutput("network",height = "750px"))
    ),
    fluidRow(
      column(12,
             tabsetPanel(
               id = 'table',
               tabPanel("nodes", DT::dataTableOutput("nodes_table")),
               tabPanel("modules", DT::dataTableOutput("module_table"))
               )
      )
      
      #column(4,
      #        h4("Node Info"), 
      #        downloadButton(outputId = "save_node_info", label = "Save Node"),
      #        htmlOutput(outputId = "node_info")
      #),
      #column(4, 
      #       h4("Module Info"),
      #       fluidRow(
      #        column(1,downloadButton(outputId ="save_module_info", label = "Save Module")),
      #        column(1, offset = 2,
      #        conditionalPanel(
      #          condition = "input.method != 'module'",
      #          downloadButton(outputId = "save_all_module_info", label = "Save all Modules"))
      #        )
      #       ),
      #       htmlOutput(outputId = "module_info")
      #       )
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
        sub_net(diffScoreSubgraph(Net, percentile()/100))
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
    if(length(gene_list()> 0 )){
    sub_net(geneListSubgraph(Net, Module, gene_list(), input$search_additional))
    }
  })
  
  
  observeEvent(input$steiner, {
    st <- buildSteinerTrees(Net, sub_net() %N>%  pull(feature))
    st <- st %E>%
      mutate(is_steiner = TRUE) %N>%
      mutate(is_steiner = TRUE)
    sub_net(graph_join(sub_net(), st) %>%
      mutate(is_steiner = replace_na(is_steiner, FALSE)) %>%
      mutate(color_code = if_else(is_steiner, "#FF0000", "#666")))
    gene_list(sub_net() %N>% pull(feature))
  })
  
  observeEvent(input$method, {
  sub_net(tbl_graph())
  })
  
  output$save_net_image <- downloadHandler(
         filename = function(){ifelse(str_length(input$file_name) > 0, paste(input$file_name, '.html', sep = ''), 'file.html')},
         content = function(file) {
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
             op <- .25
             fs <- 12
           }
           
           if(input$method == "diff"){
             S_nodes <- S_nodes %>% mutate(size = rescale(score, to = c(4, 16)))
             g <- forceNetwork(Links = S_edges, Nodes = S_nodes,
                               Source = "from", Target = "to",
                               Value = "weight", NodeID = names,
                               Group = input$group_by, Nodesize = "size", opacity = 1, opacityNoHover = op,
                               zoom = TRUE, fontSize=fs, radiusCalculation = JS("d.nodesize"),
                               charge = -10) 
           }
           else{
             g <- forceNetwork(Links = S_edges, Nodes = S_nodes,
                               Source = "from", Target = "to",
                               Value = "weight", NodeID = names,
                               Group = input$group_by, linkColour=S_edges$color_code,
                               opacity = 1, opacityNoHover = op, zoom = TRUE, fontSize=fs,
                               charge = -20)
           }
           saveWidget(g, file)
       }})
  

  
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
    S_nodes <- prepNodeTable(S_tables[[1]])
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
      op <- .25
      fs <- 12
    }
    
    if(input$method == "diff"){
      S_nodes <- S_nodes %>% mutate(size = rescale(score, to = c(4, 16)))
      forceNetwork(Links = S_edges, Nodes = S_nodes,
                   Source = "from", Target = "to",
                   Value = "weight", NodeID = names,
                   Group = input$group_by, Nodesize = "size", opacity = 1, opacityNoHover = op,
                   zoom = TRUE, fontSize=fs, radiusCalculation = JS("d.nodesize"),
                   charge = -10) 
    }
    else{forceNetwork(Links = S_edges, Nodes = S_nodes,
                 Source = "from", Target = "to",
                 Value = "weight", NodeID = names,
                 Group = input$group_by, linkColour=S_edges$color_code,
                 opacity = 1, opacityNoHover = op, zoom = TRUE, fontSize=fs,
                 charge = -20)
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

