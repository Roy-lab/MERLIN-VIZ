library(shiny)
library(tools)
library(tidyverse)
library(networkD3)

## Load in aux functions and data
load('net_data.Rdata')
source('asp_aux_functions.R')
MyClickScript <- 'Shiny.setInputValue("node_name", d.name)'


if(interactive()){
################################# ui ###########################################
    ui <- fluidPage(
    title = "AsperNet",
    
    forceNetworkOutput("network"),
    
    fluidRow(
      column(4,
             h4("Search Network"),
             selectInput(inputId = "method", "Search Method",
                  c(Modules = "module", 'GO Term' = "go_term", 'Gene' = 'gene', 'Gene List' = 'list')
             ),
             # Only show this panel if the plot type is a histogram
              conditionalPanel(
              condition = "input.method == 'module'",
                selectInput(
                inputId = "module_id", "Module ID",
                sort(module_ids)
                ),
              ),
              conditionalPanel(
              condition = "input.method == 'go_term'",
              selectInput(
                inputId = "go_term", "GO terms",
                sort(enriched_go_terms),
              ),
              radioButtons("stiener", h3("Show Stiener Tree"),
                           choices = list("Do Not Show" = 1, "Show Stiener Tree" = 2))
              ),
              conditionalPanel(
                condition = "input.method == 'list'",
                fileInput(inputId = "cell_list_file", "Cell List File"),
                conditionalPanel(
                  condition = "output.fileUploaded",
                  radioButtons("stiener", h2("Show stiener Tree"),
                               choices = list("Do Not Show" = 1, "Show Stiener Tree" = 2))
                )
              ),
             conditionalPanel(
               condition ="input.method =='gene'",
               selectizeInput(inputId = "gene", label = "Gene", choices = NULL)
               )
            ),
      column(4,
              h4("Node Info"), 
              downloadButton(outputId = "save_node_info", label = "Save Node"),
              htmlOutput(outputId = "node_info")
      ),
      column(4, 
             h4("Module Info"),
             fluidRow(
              column(1,downloadButton(outputId ="save_module_info", label = "Save Module")),
              column(1, offset = 2,
              conditionalPanel(
                condition = "input.method != 'module'",
                downloadButton(outputId = "save_all_module_info", label = "Save all Modules"))
              )
             ),
             htmlOutput(outputId = "module_info")
             )
  )
)
##################### Server Functions ########################################
server <- function(input, output, session) {
  updateSelectizeInput(session, 'gene', choices = genes, server = TRUE)
  #### Variables ################
  node_name_info <- reactiveVal(value = NA)  
  module_id_info <- reactiveVal(value = NA)
  gene_name <-reactiveVal()
  
  ##### File IO ################
  file_path <- reactive({
    file <- input$cell_list_file
    req(file)
    file <- file$datapath 
  })
  
  output$fileUploaded <- reactive({
    return(!is.null(sub_net()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  sub_net <- reactive({
    if(input$method=="module"){
      module <- input$module_id
      print(module)
      module_id_info(input$module_id)
      node_name_info(NA)
      moduleSubgraph(Module, module)
    }else if(input$method =="go_term"){
      node_name_info(NA)
      module_id_info(NA)
      goSubgraph(Net, Module, enrich_2_module, input$go_term)
    }else if(input$method =="gene"){
      geneSubgraph(Net, Module, list(input$gene))
    }
    else{
      node_name_info(NA)
      module_id_info(NA)
      fp <- file_path()
      gene_list <- read_csv(file = fp, col_names=FALSE) %>% 
        pull(X1)
      geneListSubgraph(Net, Module, gene_list)
    }
  })
  
  ############## Main Render ######################
  output$network <- renderForceNetwork({
    S<-sub_net()
    print(input$stiener)
    if(input$stiener == 1){
      S_tables <- graph2NodeEdgeTables(S)
      S_nodes <- S_tables[[1]]
      S_edges <- S_tables[[2]]
    ### Add superficial edge to overcome empty array 
    S_edges <- S_edges %>% add_row(from = 0, to = 0, weight = 0)
    
    forceNetwork(Links = S_edges, Nodes = S_nodes,
                 Source = "from", Target = "to",
                 Value = "weight", NodeID = "feature",
                 Group = "module", opacity = 1, zoom = TRUE, fontSize=40,
                 charge = -10, clickAction = MyClickScript) 
    }else{
      showModal(modalDialog(
        title = "Building Stiener Tree",
        "Stiener tree construction is underway.",
        easyClose = TRUE))
      
      Stiener <- buildStienerTrees(Net, S %N>% as_tibble() %>% pull(feature))
      
      Stiener <- Stiener %E>%
                  mutate(is_stiener = TRUE) %N>%
                  mutate(is_stiener = TRUE)
      sub_net <- graph_join(S, Stiener)
      stiener_tables <- graph2NodeEdgeTables(sub_net)
      stiener_nodes <- stiener_tables[[1]]
      stiener_edges <- stiener_tables[[2]]
      
      stiener_edges <- stiener_edges %>% 
        mutate(is_stiener = replace_na(is_stiener, FALSE)) %>%
        mutate(color_code = if_else(is_stiener, "#FF0000", "#666")) %>%
        add_row(from = 0, to = 0, weight = 0)
      
      forceNetwork(Links = stiener_edges, Nodes = stiener_nodes,
                   Source = "from", Target = "to",
                   Value = "weight", NodeID = "feature",
                   Group = "module", linkColour = stiener_edges$color_code,
                   opacity = 1, zoom = TRUE, fontSize=40,
                   charge = -10, clickAction = MyClickScript) 
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
        paste(file_path_sans_ext(basename(fp)),"_modules.txt", sep='')  
      }
    },
    content = function(con) {
      text_info <- printAllModuleInfo(sub_net(), Module, list())
      text_info <- str_replace_all(text_info, '<br/>', '\n')
      text_info <- str_replace_all(text_info, '&emsp', '\t')
      write(text_info, con)
    }
  )
  
  ############### observe Events ##########################################
  observeEvent(input$node_name, {
    node_name_info(input$node_name)
    if(input$method!="module"){
      module_id_info(getModuleID(Net, node_name_info()))
    }
  })
  
  observeEvent(node_name_info, {
    output$node_info <- renderUI({
      HTML(printNodeInfo(Net, node_name_info()))
    })
  })
  
  observeEvent(module_id_info,{
    output$module_info<-renderUI({
      HTML(printModuleInfo(Module, module_id_info(), list()))
    })
  })
  
}
########################################################################

shinyApp(ui = ui, server = server)
}
