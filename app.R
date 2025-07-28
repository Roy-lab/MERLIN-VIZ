#spencers new version that im changing but now adding new tabs

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
library(shinyWidgets)
library(shinythemes)#CC: themes package
library(shinyBS)#CC: Widgets package
library(htmltools)
library(bsplus)#CC: For tooltip
library(shinyjs)#CC: For tooltip
library(patchwork)


## Load in aux functions and data
source('aux_functions.R') 
source('printerFunction.R')
source('heatmap_printer.R')


#MyClickScript <- 'Shiny.setInputValue("save_module", module_name)'

all_gene_names <- unique(c(genes, genename_map$common_name))
palettes_nodes<- tibble(rownames_to_column(brewer.pal.info, var = 'pal')) %>% 
  filter(category == "qual")
palettes_edges <- tibble(rownames_to_column(brewer.pal.info, var = 'pal')) %>% 
  filter(category == "div")
igraph_layout <- c('Fruchterman-Reingold'='nicely', 'Davidson-Harel'='dh', 'Kamada-Kawai'='kk', 'Large graph layout'= 'lgl') #'Force directed' = 'drl')

################################# ui ###########################################
ui <- navbarPage(title,
                 id = 'navbar',
                 theme = shinytheme("flatly"), #CC: Set this as the main theme
                 
                 tabPanel("Visualize",
                          #shinythemes::themeSelector(), #CC: this allows you to cycle through different themes when the app is running
                          #####  Gene selection Options ##############################
                          fluidRow(
                            column(2,
                                   pickerInput(inputId = "method", h4("Search Method"), #CC: I set the default to list
                                               c('Gene List' = 'list', Modules = "module", "Node Diffusion" = "diff",'GO-Term' = "go_term"),
                                               selected = "list"
                                   ),
                                   #### Module Selection ###################
                                   conditionalPanel(
                                     condition = "input.method == 'module'",
                                     pickerInput(
                                       inputId = "module_id", 
                                       label = "Module ID",
                                       choices = c("", unlist(sort(module_ids)))
                                     )
                                   ),
                                   #### GO Selection  ###################
                                   conditionalPanel(
                                     condition = "input.method == 'go_term'",
                                     pickerInput(inputId = "go_term", label = "GO terms", choices = c("", unlist(sort(enriched_go_terms))))
                                   ),
                                   
                                   ### Gene List Selection ###############
                                   conditionalPanel(
                                     condition = "input.method == 'list'",
                                     tags$div(
                                       style = "display: flex; align-items: center;",
                                       tags$h4("Input Genes"),
                                       tags$div(
                                         style = "margin-left: 1px;", # CC: moved make the icon and the label closer
                                         bsButton("Inputgenes", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), #CC: tooltip for Gene selection. Note that itll stop working if the sentence is too long
                                         bsPopover("Inputgenes", "Additional Info",
                                                   "Select your genes of interest. Those ending with _NCA represent transcription factors whose activities are based on binding motifs, rather than gene expression",
                                                   placement = "right",
                                                   options = list(
                                                     container = "body",
                                                     html = TRUE,
                                                     template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                   )
                                         )
                                       )
                                     ),
                                     selectInput(inputId = "gl", label = NULL, choices = NULL, multiple = TRUE, selectize = TRUE, selected = c("srbA")),
                                     
                                     tags$div(
                                       style = "display: flex; align-items: center;",
                                       tags$h4("Or Upload Gene List"),
                                       tags$div(
                                         style = "margin-left: 1px;",
                                         bsButton("uList", "", icon = icon("question-circle", class = "fa-lg"), style = "link"),
                                         bsPopover("uList", "Additional Info",
                                                   "Upload text file with AFUA gene names. One gene per row",
                                                   placement = "right",
                                                   options = list(
                                                     container = "body",
                                                     html = TRUE,
                                                     template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                   )
                                         )
                                       )
                                     ),
                                     fileInput(inputId = "cell_list_file", "")
                                   ),
                                   
                                   ############# Gene List Selection  2 ###########################
                                   conditionalPanel(
                                     condition = "input.method == 'list'",
                                     checkboxGroupInput(
                                       inputId = "search_additional",
                                       label = tags$div(
                                         style = "display: flex; align-items: center;",
                                         tags$h4("Additional Options"),
                                         tags$div(
                                           style = "margin-left: 1px;", 
                                           bsButton("Seachinfo", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                           bsPopover("Seachinfo", "Additional Info",
                                                     "Neighbors are genes that have a direct connection to your query; Module members share the same regulatory program. A Steiner tree finds the smallest path between two genes.",
                                                     placement = "right",
                                                     options = list(
                                                       container = "body",
                                                       html = TRUE,
                                                       template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                     )
                                           )
                                         )
                                       ),
                                       choices = c("Show Neighbors" = "neigh", "Include Module Members" = "mod", "Create Steiner Tree" = "stein"), #CC: Changed to more descriptive titles
                                       selected = c("mod")
                                      )
                                   ),
                                   #This is required to align the tooltip and the download button in the Diffusion Score File section
                                   tags$style(HTML("
   /* Custom CSS for the 'Example' button */
    #download_example_diff {
     padding: 5px 7px;
      font-size: 9px;   
      margin-left: -50px;
      margin-top: 5px;
    }
    /* Custom CSS for horizontal alignment */
    .horizontal-align {
      display: flex;
      align-items: center;
    
    }
    

  ")),
                                   
                                   ####### Diffusion Selection ##################
                                   conditionalPanel(
                                     condition = "input.method == 'diff'",
                                     fluidRow(
                                       column(
                                         width = 9,
                                         tags$h4("Score File"),
                                         fileInput(inputId = "diff_list_file", ""),
                                         verbatimTextOutput("example_text") #CC: added a an example file for diffuion analysis. 
                                       ),
                                       column(
                                         width = 2,
                                         style = "margin-top: -8px; margin-left: 1px;",
                                         div(
                                           class = "horizontal-align",  # Apply the CSS class for horizontal alignment
                                           downloadButton("download_example_diff",""),
                                           bsButton("additional_info_diff", "", icon = icon("question-circle", class = "fa-lg"), style = "link")
                                         ),
                                         bsPopover("additional_info_diff", "Additional Info",
                                                   "A diffuision analysis requires Genes in the 1st column, and any numerical value tied to that gene in the 2nd column; such as the log P-value.",
                                                   placement = "right",
                                                   options = list(
                                                     container = "body",
                                                     html = TRUE,
                                                     template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                   )
                                         ),
                                         bsPopover("download_example_diff", "Additional Info",
                                                   "An example file containg differentially expressed genes in the 1st column, and the log P-value of each gene in the 2nd column. From Rush et al., 2019.",
                                                   placement = "right",
                                                   options = list(
                                                     container = "body",
                                                     html = TRUE,
                                                     template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                   )
                                         )
                                       )
                                     ),
                                     
                                     div(
                                       numericInput(
                                         inputId = "min_neigh",
                                         label = tags$div(
                                           style = "display: flex; align-items: center;",
                                           tags$h4("Min # of Targets"),
                                           tags$div(
                                             style = "margin-left: 1px;",
                                             bsButton("additional_info_min_neigh", "", icon = icon("question-circle", class = "fa-lg"), style = "link"),
                                             bsPopover("additional_info_min_neigh", "Additional Info",
                                                       "Specify the minimum number of target genes for diffusion analysis.",
                                                       placement = "right",
                                                       options = list(
                                                         container = "body",
                                                         html = TRUE,
                                                         template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                       )
                                             )
                                           )
                                         ),
                                         value = 5
                                       )
                                     ),
                                     
                                     div(
                                       selectInput(
                                         inputId = "kernel",
                                         label = tags$div(
                                           style = "display: flex; align-items: center;",
                                           tags$h4(" Lambda score"),
                                           tags$div(
                                             style = "margin-left: 1px;",
                                             bsButton("Lambda", "", icon = icon("question-circle", class = "fa-lg"), style = "link"),
                                             bsPopover("Lambda", "Additional Info",
                                                       "Lambda refers to the laplacian kernel diffusion constant. Larger lambda increase diffusion distance resulting in a smoother resulting score.",
                                                       placement = "right",
                                                       options = list(
                                                         container = "body",
                                                         html = TRUE,
                                                         template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                       )
                                             )
                                           )
                                         ),
                                         c(1, 10, 100, 1000)
                                       )
                                     ),
                                     ###
                                     div(
                                       numericInput(
                                         inputId = "disp_regs",
                                         label = tags$div(
                                           style = "display: flex; align-items: center;",
                                           tags$h4("# of Regulators to Display"),
                                           tags$div(
                                             style = "margin-left: 1px;",
                                             bsButton("Number of Regulators to Display", "", icon = icon("question-circle", class = "fa-lg"), style = "link"),
                                             bsPopover("Number of Regulators to Display", "Additional Info",
                                                       "Here you can specify the number of regulators for diffusion analysis.",
                                                       placement = "right",
                                                       options = list(
                                                         container = "body",
                                                         html = TRUE,
                                                         template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                       )
                                             )
                                           )
                                         ), 
                                         value = 5
                                       )
                                     ),
                                     
                                     actionButton(inputId = "refresh_diff", "Refresh")
                                     
                                   ),
                                   # br(),
                                   # actionButton("contact_button", "Contact us",   
                                   # tags$script(HTML("
                                   #      // JavaScript code to open a link
                                   #         document.getElementById('go_to_page').onclick = function() {
                                   #           window.open('https://www.example.com', '_blank');
                                   #          }")))
                            ),
                            #### Visualization Block ########
                            column(10,
                                   tabsetPanel(id = "displayType", type = "tabs",
                                               tabPanel("Network Plot", plotOutput("print_net", click = 'plot_click', height = '1000px')), #CC: Changed to more descriptive titles
                                               tabPanel("Expression Heatmaps", plotOutput("expression_heatmap", height = '1000px')),
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
                                     selectInput(inputId = 'print_layout', choices = igraph_layout, multiple = FALSE, 
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Node layout"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("NodeLayout", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("NodeLayout", "Additional Info",
                                     	    			  "Default layouts for node display. See igraph package for more details of each method.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     ))),
                                     sliderInput(inputId = 'print_min_genes', value = 1, min = 1, max = 10,
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Minimum number of genes in component"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("CCComps", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("CCComps", "Additional Info",
                                     	    			  "The minimum number of genes to be contained in a connected component to display. e.g If this is set to 2, then genes that are have no neighbors will be removed from display.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     ))),
                                     selectInput(inputId = 'print_disp_names', choices = NULL, multiple = TRUE, selectize = TRUE,
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Display gene names"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("DispNames", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("DispNames", "Additional Info",
                                     	    			  "The set of genes to display with labels. Toggle a gene either by clicking on it in the display or by addition to this.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     	    	))),
                                     radioButtons("common_name",
                                                  choices = list("Common" = 1, "Systematic" = 2),
                                     		label = tags$div(
                                     	     	style = "display: flex; align-items: center;",
                                     	     	tags$h4("Name format"),
                                     	     	tags$div(
                                     	     		style = "margin-left: 1px;", 
                                     	     		bsButton("NFInfo", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	     		bsPopover("NFInfo", "Additional Info",
                                     	     			  "Name format for display. Select common names to use first instance of name from fungiDB database. Systematic names are in AFUA_#G##### format.",
                                     	     			  placement = "right",
                                     	     			  options = list(
                                     	     			  	container = "body",
                                     	     			  	html = TRUE,
                                     	     			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	     			  )
                                     	     		)
                                     	     	))),
                                     checkboxInput(inputId = 'print_name_bool', label = "Name in node", value = TRUE),
                                     conditionalPanel(
                                       condition = "input.print_name_bool == 0",
                                       sliderInput(inputId = 'print_nudge_y', value = 0.3, min = 0, step = 0.1, max = 5, 
                                       	    label = tags$div(
                                       	    	style = "display: flex; align-items: center;",
                                       	    	tags$h4("Nudge labels"),
                                       	    	tags$div(
                                       	    		style = "margin-left: 1px;", 
                                       	    		bsButton("y_nudge_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                       	    		bsPopover("y_nudge_info", "Additional Info",
                                       	    			  "Move node labels (y axis).",
                                       	    			  placement = "right",
                                       	    			  options = list(
                                       	    			  	container = "body",
                                       	    			  	html = TRUE,
                                       	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                       	    			  )
                                       	    		)
                                       	    	))),
                                       sliderInput(inputId = 'print_text_angle', value = 0, min = -90, max = 90, step = 15, 
                                       	    label = tags$div(
                                       	    	style = "display: flex; align-items: center;",
                                       	    	tags$h4("Text angle"),
                                       	    	tags$div(
                                       	    		style = "margin-left: 1px;", 
                                       	    		bsButton("name_angle_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                       	    		bsPopover("name_angle_info", "Additional Info",
                                       	    			  "Rotate node labels.",
                                       	    			  placement = "right",
                                       	    			  options = list(
                                       	    			  	container = "body",
                                       	    			  	html = TRUE,
                                       	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                       	    			  )
                                       	    		)
                                       	    	)))
                                     )
                              ),
                              column(2,
                                     radioButtons(inputId = "print_group_by",
                                                  choices = list ("Module" = 'module', "Regulator" = "regulator", "Gene Name" = "geneSuper"), selected =  'module', 
                                     	     label = tags$div(
                                     	     	style = "display: flex; align-items: center;",
                                     	     	tags$h4("Node color by"),
                                     	     	tags$div(
                                     	     		style = "margin-left: 1px;", 
                                     	     		bsButton("color_by_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	     		bsPopover("color_by_info", "Additional Info",
                                     	     			  "Node coloring method. If set to module, nodes are colored by module assigment. Grey nodes correspond to genes not assigned to a module. If set to regulator, regulators are colored red and targets are colored blue. If set to Gene Name, color nodes with similar common gene name as the same color.",
                                     	     			  placement = "right",
                                     	     			  options = list(
                                     	     			  	container = "body",
                                     	     			  	html = TRUE,
                                     	     			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	     			  )
                                     	     		)
                                     	     	))), #CC: Removed gene family labeling as it was confusing to users
                                     selectInput(inputId = 'print_node_pal', choices = palettes_nodes$pal, multiple = FALSE, selected = 'Pastel1', 
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Node color palette"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("pal_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("pal_info", "Additional Info",
                                     	    			  "Palette options to color nodes in display field. Options are provided by the color brewer package.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     	    	))),
                                     sliderInput(inputId = 'print_max_node_size', value = 8, min = 1, max = 25, 
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Node size"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("node_size_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("node_size_info", "Additional Info",
                                     	    			  "The display size of nodes that do not contain text. Nodes that contain text will be scaled to compensate for text size.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     	    	))),
                                     sliderInput(inputId = 'print_font_size', value = 8, min = 1, max = 25, 
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Node label font size"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("node_text_size_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("node_text_size_info", "Additional Info",
                                     	    			  "The font size of node labels.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     	    	)))
                                     
                              ),
                              column(2,
                                     radioButtons(inputId = "edge_color_by", 
                                                  choices = list("Correlation" = "Correlation", "Regression Weight"= "Reg_weight"), label = tags$div(
                                                    style = "display: flex; align-items: center;",
                                                    tags$h4("Edge color by"),
                                                    tags$div(
                                                      style = "margin-left: 1px;", 
                                                      bsButton("edge_color_by_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                      bsPopover("edge_color_by_info", "Additional Info",
                                                                "Edge coloring method. Edges will be colored by either regression weight or correlation betweeen linked genes.",
                                                                placement = "right",
                                                                options = list(
                                                                  container = "body",
                                                                  html = TRUE,
                                                                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                                )
                                                      )
                                                    ))),
                                     
                                     conditionalPanel(
                                       condition = "input.edge_color_by == 'Reg_weight'",
                                       sliderInput("edge_color_range", label = tags$div(
                                         style = "display: flex; align-items: center;",
                                         tags$h4("Edge color range"),
                                         tags$div(
                                           style = "margin-left: 1px;", 
                                           bsButton("edge_color_range_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                           bsPopover("edge_color_range_info", "Additional Info",
                                                     "Set the minimum and maximum of the regression weight color scale.",
                                                     placement = "right",
                                                     options = list(
                                                       container = "body",
                                                       html = TRUE,
                                                       template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                     )
                                           ))), min = -10, max = 10, value = c(-5,5), step = 0.1),
                                     ),
                                     selectInput(inputId = 'edge_color_palette', choices = palettes_edges$pal, multiple = FALSE, selected = 'RdBu', 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4("Edge color palette"),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("edge_pal_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("edge_pal_info", "Additional Info",
                                                               "Palette options to color edges in display field. Options are provided by the color brewer package.",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   ))),
                                     ),
                              column(2,
                                     sliderInput(inputId = 'print_expand_x', value = 2, step = 0.1, min = 0, max = 25, 
                                     	    	    label = tags$div(
                                     	    	    	style = "display: flex; align-items: center;",
                                     	    	    	tags$h4("Expand X axis"),
                                     	    	    	tags$div(
                                     	    	    		style = "margin-left: 1px;", 
                                     	    	    		bsButton("expand_X_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    	    		bsPopover("expand_X_info", "Additional Info",
                                     	    	    			  "Expands X axis. Scale this to fit large node names in plot window.",
                                     	    	    			  placement = "right",
                                     	    	    			  options = list(
                                     	    	    			  	container = "body",
                                     	    	    			  	html = TRUE,
                                     	    	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    	    			  )
                                     	    	    		)
                                     	    	    	))),
                                     sliderInput(inputId = 'print_expand_y', value = 0, step = 0.1, min = 0, max = 25, 
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Expand Y axis"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("expand_Y_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("expand_Y_info", "Additional Info",
                                     	    			  "Expands Y axis. Scale this to fit large node names in plot window.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     	    		))),
                                     	    		sliderInput("legend_font_size", label = tags$div(
                                     	    		  style = "display: flex; align-items: center;",
                                     	    		  tags$h4("Legend font size"),
                                     	    		  tags$div(
                                     	    		    style = "margin-left: 1px;", 
                                     	    		    bsButton("legend_font_size_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		    bsPopover("legend_font_size_info", "Additional Info",
                                     	    		              "Set the legend font size of the figure. Note if too small, the font will not display.",
                                     	    		              placement = "right",
                                     	    		              options = list(
                                     	    		                container = "body",
                                     	    		                html = TRUE,
                                     	    		                template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    		              )
                                     	    		    ))), min = 0, max = 25, value = 18, step = 1),
                                     sliderInput(inputId = 'print_image_height', value = 8, min = 1, max = 25, 
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4("Image height (in)"),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("image_height_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("image_height_info", "Additional Info",
                                     	    			  "Set height (in inches) of image when saved. Save by hitting the save figure button.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     	    	))),
                                     sliderInput(inputId = 'print_image_width', value = 8, min = 1, max = 25, 
                                     	    label = tags$div(
                                     	    	style = "display: flex; align-items: center;",
                                     	    	tags$h4('Image width (in)'),
                                     	    	tags$div(
                                     	    		style = "margin-left: 1px;", 
                                     	    		bsButton("image_width_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     	    		bsPopover("image_width_info", "Additional Info",
                                     	    			  "Set width (in inches) of image when saved. Save by hitting the save figure button.",
                                     	    			  placement = "right",
                                     	    			  options = list(
                                     	    			  	container = "body",
                                     	    			  	html = TRUE,
                                     	    			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     	    			  )
                                     	    		)
                                     	    	))),
                                     textInput("print_file_name", value = "file_name", 
                                     	label = tags$div(
                                     	style = "display: flex; align-items: center;",
                                     	tags$h4('File name'),
                                     	tags$div(
                                     		style = "margin-left: 1px;", 
                                     		bsButton("figure_name_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                     		bsPopover("figure_name_info", "Additional Info",
                                     			  "Set custom figure name for when figure is saved.",
                                     			  placement = "right",
                                     			  options = list(
                                     			  	container = "body",
                                     			  	html = TRUE,
                                     			  	template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                     			  )
                                     		)
                                     	))
                                     	  ),
                                     downloadButton("saveFig", "Save figure")
                              )
                            ),
                            ############# Heatmap Suport option Tab  ########################
                            conditionalPanel(
                              condition ="input.displayType == 'Expression Heatmaps'",
                              column(2),
                              column(2,
                                     radioButtons("common_name_heatmap",
                                                  choices = list("Common" = 1, "Systematic" = 2),
                                                  label = tags$div(
                                                    style = "display: flex; align-items: center;",
                                                    tags$h4("Name format"),
                                                    tags$div(
                                                      style = "margin-left: 1px;", 
                                                      bsButton("NFInfo2", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                      bsPopover("NFInfo2", "Additional Info",
                                                                "Name format for display. Select common names to use first instance of name from fungiDB database. Systematic names are in AFUA_#G##### format.",
                                                                placement = "right",
                                                                options = list(
                                                                  container = "body",
                                                                  html = TRUE,
                                                                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                                )
                                                      )
                                                    )))
                              ),
                              column(2,
                                     selectInput(inputId = 'tfa_palette_heatmap', choices = palettes_edges$pal, multiple = FALSE, selected = 'PiYG', 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4("TFA color palette"),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("expr_pal_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("expr_pal_info", "Additional Info",
                                                               "Palette to color transcription factor activity profile heatmaps.",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   ))),
                                     sliderInput(inputId = 'tfa_range_heatmap', value = c(-2,2), min = -10, max = 10, 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4("TFA range"),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("tfa_range_heatmap_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("tfa_range_heatmap_info", "Additional Info",
                                                               "Range to color TFA profile heatmap",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   ))),
                                     selectInput(inputId = 'expression_palette_heatmap', choices = palettes_edges$pal, multiple = FALSE, selected = 'RdBu', 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4("Expression color palette"),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("expression_palette_heatmap_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("expression_palette_heatmap_info", "Additional Info",
                                                               "Palette options to color expression heatmap",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   ))),
                                     sliderInput(inputId = 'expression_range_heatmap', value = c(-2,2), min = -10, max = 10, 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4("Expression range"),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("expression_range_heatmap_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("expression_range_heatmap_info", "Additional Info",
                                                               "Range to color expression heatmap",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   )))
                                     
                              ),
                              column(2,
                                     radioButtons(inputId = "edge_color_by_heatmap", 
                                                  choices = list("Correlation" = "Correlation", "Regression Weight"= "Reg_weight"), label = tags$div(
                                                    style = "display: flex; align-items: center;",
                                                    tags$h4("Edge color by"),
                                                    tags$div(
                                                      style = "margin-left: 1px;", 
                                                      bsButton("edge_color_by_info2", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                      bsPopover("edge_color_by_info2", "Additional Info",
                                                                "Edge coloring method. Edges will be colored by either regression weight or correlation betweeen linked genes.",
                                                                placement = "right",
                                                                options = list(
                                                                  container = "body",
                                                                  html = TRUE,
                                                                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                                )
                                                      )
                                                    ))),
                                     conditionalPanel(
                                       condition = "input.edge_color_by_heatmap == 'Reg_weight'",
                                       sliderInput("edge_color_range_heatmap", label = tags$div(
                                         style = "display: flex; align-items: center;",
                                         tags$h4("Edge color range"),
                                         tags$div(
                                           style = "margin-left: 1px;", 
                                           bsButton("edge_color_range_info2", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                           bsPopover("edge_color_range_info2", "Additional Info",
                                                     "Set the minimum and maximum of the regression weight color scale.",
                                                     placement = "right",
                                                     options = list(
                                                       container = "body",
                                                       html = TRUE,
                                                       template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                     )
                                           ))), min = -10, max = 10, value = c(-5,5), step = 0.1),
                                     ),
                                     selectInput(inputId = 'edge_color_palette_heatmap', choices = palettes_edges$pal, multiple = FALSE, selected = 'RdBu', 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4("Edge color palette"),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("edge_pal_info2", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("edge_pal_info2", "Additional Info",
                                                               "Palette options to color edges in display field. Options are provided by the color brewer package.",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   )))
                                    
                              ),
                              column(2,
                                     sliderInput("Font_size_heatmap", label = tags$div(
                                       style = "display: flex; align-items: center;",
                                       tags$h4("Figure font size"),
                                       tags$div(
                                         style = "margin-left: 1px;", 
                                         bsButton("legend_font_size_info2", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                         bsPopover("legend_font_size_info2", "Additional Info",
                                                   "Set the legendfont size of the figure. Note if too small, the font will not display.",
                                                   placement = "right",
                                                   options = list(
                                                     container = "body",
                                                     html = TRUE,
                                                     template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                   )
                                         ))), min = 0, max = 25, value = 18, step = 1),
                                     sliderInput(inputId = 'print_image_height_heatmap', value = 8, min = 1, max = 25, 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4("Image height (in)"),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("image_height_info2", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("image_height_info2", "Additional Info",
                                                               "Set height (in inches) of image when saved. Save by hitting the save figure button.",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   ))),
                                     sliderInput(inputId = 'print_image_width_heatmap', value = 8, min = 1, max = 25, 
                                                 label = tags$div(
                                                   style = "display: flex; align-items: center;",
                                                   tags$h4('Image width (in)'),
                                                   tags$div(
                                                     style = "margin-left: 1px;", 
                                                     bsButton("image_width_info2", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                     bsPopover("image_width_info2", "Additional Info",
                                                               "Set width (in inches) of image when saved. Save by hitting the save figure button.",
                                                               placement = "right",
                                                               options = list(
                                                                 container = "body",
                                                                 html = TRUE,
                                                                 template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                               )
                                                     )
                                                   ))),
                                     textInput("print_file_name_heatmap", value = "file_name", 
                                               label = tags$div(
                                                 style = "display: flex; align-items: center;",
                                                 tags$h4('File name'),
                                                 tags$div(
                                                   style = "margin-left: 1px;", 
                                                   bsButton("figure_name_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                                                   bsPopover("figure_name_info", "Additional Info",
                                                             "Set custom figure name for when figure is saved.",
                                                             placement = "right",
                                                             options = list(
                                                               container = "body",
                                                               html = TRUE,
                                                               template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                                                             )
                                                   )
                                                 ))
                                     ),
                                     downloadButton("saveFig2", "Save figure")
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
          pal <- extend_palette <- colorRampPalette(pal)(num_mods)
        }
        colorScale = JS(paste0('color=d3.scaleOrdinal([ `#BEBEBE`, ', paste(sprintf('`%s`', pal), collapse = ', '), ']), color.domain([-9999])'))
      }else if(input$group_by == "geneSuper"){
        num_supers <- length(setdiff(unique(S_nodes %>% pull(geneSuper)), "Unlabeled"))
        max_colors <- palettes$maxcolors[which(palettes$pal == "Pastel1")]
        pal <- brewer.pal(max(3, min(max_colors, num_supers)), "Pastel1")
        if(max_colors < num_supers){
          pal <- extend_palette <- colorRampPalette(pal)(num_supers)
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
        makeSubNetGraph(subNet, names_in_nodes = input$print_name_bool, node_color_by = input$print_group_by, 
                        edge_color_by = input$edge_color_by, edge_color_palette = input$edge_color_palette, 
                        node_color_palette = input$print_node_pal, 
                        node_size_by = node_size_by, max_node_size = input$print_max_node_size, 
                        layout = input$print_layout, focus_nodes = list(), 
                        font_size = input$print_font_size, 
                        nudge_y = input$print_nudge_y, text_angle = input$print_text_angle, show_legend = TRUE,
                        expand_x = input$print_expand_x, expand_y = input$print_expand_y, color_scale_limits = input$edge_color_range, legend_font_size = input$legend_font_size)
      )
      gg_out_plot()
    }
  })
  
  
  output$expression_heatmap <- renderPlot({
    subNet <- sub_net() 
    if(is_empty(subNet)){
      gg<- ggplot() + 
        theme(panel.background = element_rect(fill="white", colour = "white")) +
        geom_text(label = "no subgraph selected.")
    }else{
      print(input$tfa_color_range_heatmap)
      gg_out_heatmap(
        makeSubgraphHeatmap(subNet,
          display_name = input$common_name_heatmap,
          edge_color_by = input$edge_color_by_heatmap, 
          edge_color_palette = input$edge_color_palette_heatmap, 
          font_size = input$print_font_size,
          tfa_color_palette = input$tfa_palette_heatmap, 
          expression_color_palette = input$expression_palette_heatmap,
          scale_edge_color  = input$edge_color_range_heatmap, 
          scale_expression_colors = input$expression_range_heatmap, 
          scale_tfa_colors = input$tfa_range_heatmap, 
          figure_font_size = input$Font_size_heatmap ###_____TTEMPT######
        )
      )
      gg_out_heatmap()
    }
  })
  
  
  
  
  
  output$saveFig <- downloadHandler(
    filename = function(){ifelse(str_length(input$print_file_name) > 0, paste(input$print_file_name, '.png', sep = ''), 'file.png')},
    content = function(file){
      ggsave(file,gg_out_plot(), width = input$print_image_width, height = input$print_image_height,units = 'in')
    })
  
  
  output$saveFig2 <- downloadHandler(
    filename = function(){ifelse(str_length(input$print_file_name) > 0, paste(input$print_file_name, '.png', sep = ''), 'file.png')},
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


