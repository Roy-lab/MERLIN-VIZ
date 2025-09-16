## UI control options, these used to refactor and clean up the interface section
#of app.R to improve readability and functionality. 

# Left Banner Search IO Implementation -----

## Module Search Method -----
### Method Picker -----
method_Picker <- function(){
  pickerInput(inputId = "method", h4("Search Method"), #CC: I set the default to list
            c('Gene List' = 'list', Modules = "module", "Node Diffusion" = "diff",'GO-Term' = "go_term"),
            selected = "list"
)
}

### Module Picker -----
module_Picker <- function(module_ids){
  pickerInput(
  inputId = "module_id", 
  label = "Module ID",
  choices = c("", unlist(sort(module_ids)))
)
}



## Gene List Search Method -----
### Gene List Picker ---- 
geneList_Picker_Tag <- function(){
  tags$div(
    style = "display: flex; align-items: center;",
    tags$h4("Input Genes"),
    tags$div(
      style = "margin-left: 1px;", # CC: moved make the icon and the label closer
      bsButton("Inputgenes", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
      bsPopover("Inputgenes", "Additional Info",
                "Select your genes of interest.",
                placement = "right",
                options = list(
                  container = "body",
                  html = TRUE,
                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                )
      )
    )
  )
}


geneList_Picker <- function(){
  selectInput(inputId = "gl", label = NULL, choices = NULL, multiple = TRUE, selectize = TRUE, selected = c(""))
}


### Gene List File Input ----
geneList_File_Tag <- function(){
  tags$div(
    style = "display: flex; align-items: center;",
    tags$h4("Or Upload Gene List"),
    tags$div(
      style = "margin-left: 1px;",
      bsButton(
        "uList",
        "",
        icon = icon("question-circle", class = "fa-lg"),
        style = "link"
      ),
      bsPopover(
        "uList",
        "Additional Info",
        "Upload text file with Sytematic gene names. One gene per row",
        placement = "right",
        options = list(
          container = "body",
          html = TRUE,
          template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
        )
      )
    )
  )
}

geneList_File <- function(){
  fileInput(inputId = "cell_list_file", "")
}

### Gene List Additional Options CheckBox -----
geneList_checkBox <- function(){
  checkboxGroupInput(
    inputId = "search_additional",
    label = tags$div(
      style = "display: flex; align-items: center;",
      tags$h4("Additional Options"),
      tags$div(
        style = "margin-left: 1px;", 
        bsButton("Seachinfo", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
        bsPopover("Seachinfo", "Additional Info",
                  "Neighbors are genes that have a direct connection to your query; Module members share the same regulatory program. A Steiner tree finds the smallest (estimated) tree connecteing genes.",
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
    selected = c("neigh")
  )
}



## Diffusion Search Method -----
### Diffusion File Input ----
diffusion_File_Tag <- function(){
  tags$div(
    style = "display: flex; align-items: center;",
    tags$h4("Score File"),
    tags$div(
      style = "margin-left: 1px;",
      bsButton(
        "additional_info_diff",
        "",
        icon = icon("question-circle", class = "fa-lg"),
        style = "link"
      ),
      bsPopover("additional_info_diff", "Additional Info",
                "A diffusion analysis requires Genes in the 1st column, and any numerical value tied to that gene in the 2nd column; such as the -log P-value.",
                placement = "right",
                options = list(
                  container = "body",
                  html = TRUE,
                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                )
      )
    )
  )
}
diffusion_File <- function(){
  fileInput(inputId = "diff_list_file", "")
}



### Diffusion Min Target Input ----
diffusion_minTarget_Numerical <- function() {
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
}

### Diffusion Lambda Input ----- 
diffusion_Lambda_Select <- function() {
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
}

### Diffusion Top Regulator Input ----
diffusion_topRegulator_Numeric <- function(){
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
}

### Diffusion Refresh Action Button ---- 
diffusion_refresh_action <- function(){
  actionButton(inputId = "refresh_diff", "Refresh")
}



# Network visualization support options ------
## Graph features ------
### Print Layout Select -----
networkViz_layout_Select <- function(igraph_layout){
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
                )))
}

### Min Gene CC Slider ----- 
networkViz_minGeneCC_Slider <- function(){
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
              )))
}

### Display Name Select ----
networkViz_dispName_Select <- function(){
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
                )))
}


### Name Format Radio----
networkViz_nameFormat_Radio <- function(){
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
                 )))
}

### Name Position Bool -----
networkViz_namePositionBool_checkBox <- function(){
  checkboxInput(inputId = 'print_name_bool', label = "Name in node", value = FALSE)
}

### Name Nudge Y Slider -----
networkViz_nudgeY_Slider <- function(){
  sliderInput(inputId = 'print_nudge_y', value = 0.1, min = 0, step = 0.1, max = 5, 
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
                )))
}

### Name Angle Slider -----
networkViz_nameAngle_Slider <- function(){
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
}


## Node features -----
### Node Color Radio -----
networkViz_nodeColor_Radio <- function(){
  radioButtons(inputId = "print_group_by",
               choices = list ("Expression" = "exp", "Module" = 'module', "Regulator" = "regulator"), selected =  'exp', 
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
    )))
}

### Node Node Expression Global Checkbox -----
networkViz_nodeExpGlobal_Checkbox <- function(){
  checkboxInput("global", TRUE, label = tags$div(
    style = "display: flex; align-items: center;",
    tags$h5("All cells"),
    tags$div(
      style = "margin-left: 1px;", 
      bsButton('node_ExpGlobal', "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
      bsPopover("node_ExpGlobal", "Additional Info",
                "If selected all cells expression profiles are used to compute mean expression and edge weight. Unselected allows for coloring by cell cluster and type.",
                placement = "right",
                options = list(
                  container = "body",
                  html = TRUE,
                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                )
      ))))
}

### Node exp display cluster ----- 
networkViz_dispCluster_Select <- function(cluster){
  selectInput("disp_cluster", choices = cluster, selected = cluster[[1]], multiple = FALSE, label = tags$div(
    style = "display: flex; align-items: center;",
    tags$h5("Cluster"),
    tags$div(
      style = "margin-left: 1px;", 
      bsButton("node_Disp_Cluster", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
      bsPopover("node_Disp_Cluster", "Additional Info",
                "Cell cluster to use to compute mean expression value and edge weight.",
                placement = "right",
                options = list(
                  container = "body",
                  html = TRUE,
                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                )
      ))))
}

### Node exp disp cell type ----
networkViz_dispType_Select <- function(type){
  selectInput("disp_type", choices = type, selected = type[[1]], multiple = FALSE, label = tags$div(
  style = "display: flex; align-items: center;",
  tags$h5("Type"),
  tags$div(
    style = "margin-left: 1px;", 
    bsButton("node_Disp_Type", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
    bsPopover("node_Disp_Type", "Additional Info",
              "Cell type to compute mean expression value and edge weight.",
              placement = "right",
              options = list(
                container = "body",
                html = TRUE,
                template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
              )
    ))))
}

### Node exp color scale slider -----
networkViz_expColor_Slider <- function(){ 
sliderInput("exp_color_scale", min = 0, max = 10, value = 10, , step = 0.1, label = tags$div(
  style = "display: flex; align-items: center;",
  tags$h5("Expression Color Scale"),
  tags$div(
    style = "margin-left: 1px;", 
    bsButton("exp_color_scale_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
    bsPopover("exp_color_scale_info", "Additional Info",
              "Scale for expression used to color nodes.",
              placement = "right",
              options = list(
                container = "body",
                html = TRUE,
                template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
              )
    ))))
}

### Node color palette select ----
networkViz_nodeColor_Select<- function(palettes_nodes){
  selectInput(inputId = 'print_node_pal', choices = palettes_nodes$pal, multiple = FALSE, selected = default_node_color_pallette, 
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
                )))
}

### Node max size slider -----
networkViz_nodeSize_Slider <- function(){
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
                )))
} 

### Node font size slider -----
networkViz_nodeFontSize_Slider <- function(){
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
} 


## Edge features -----
### Edge color radio ----- 
networkViz_edgeColor_Radio <- function(){
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
                 )))
}

### Edge  color range regression slider -----
networkViz_edgeRangeReg_Slider <- function(){
  sliderInput("edge_color_range_reg", label = tags$div(
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
      ))), min = -10, max = 10, value = c(-5,5), step = 0.1)
}

### Edge color range correlation slider -----
networkViz_edgeRangeCorr_Slider <- function(){
  sliderInput("edge_color_range_corr", label = tags$div(
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
      ))), min = -1, max = 1, value = c(-1,1), step = 0.1)
}

### Edge color range abs correlation slider -----
networkViz_edgeRangeAbsCorr_Slider <- function(){
  sliderInput("edge_range_abs_corr", label = tags$div(
    style = "display: flex; align-items: center;",
    tags$h4("Edge width"),
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
      ))), min = 0, max = 4, value = 1, step = 0.1)
}

### Edge width slider ----- 
networkViz_edgeWidth_Slider <- function(){
  sliderInput("arrow_size", label = tags$div(
    style = "display: flex; align-items: center;",
    tags$h4("Arrow size"),
    tags$div(
      style = "margin-left: 1px;", 
      bsButton("arrow_size_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
      bsPopover("arrow_size_info", "Additional Info",
                "Scale arrow size ",
                placement = "right",
                options = list(
                  container = "body",
                  html = TRUE,
                  template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                )
      ))), min = 0, max = 2, value = 0.5, step = 0.1)
}

### Edge color palette select ------
networkViz_edgePalette_Select <- function(palettes_edges){
  selectInput(inputId = 'edge_color_palette', choices = palettes_edges$pal, multiple = FALSE, selected = default_edge_color_pallette, 
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
                )))
}







## Figure Save features ----
### expand X slider  ----- 
networkViz_expandX_Slider <- function(){
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
                )))
}


### expand Y slider -----
networkViz_expandY_Slider <- function(){
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
                )))
}


### Legend font size slider -----
networkViz_legendSize_Slider <- function(){
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
      ))), min = 0, max = 25, value = 18, step = 1)
}


### Image save height slider ------
networkViz_imageSaveHeight_Slider <- function(){
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
                )))

}

### Image save width slider -----
networkViz_imageSaveWidth_Slider <- function(){
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
                )))
  
}

### Image file name text -----
networkViz_imageName_Text <- function(){
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
                )))
  
}

### Image file type select -------
networkViz_imageType_Select <- function(){
  selectInput("print_file_type", "File format:",
              choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
              selected = "png")
  
}

### Save figure download -------
networkViz_Download <- function(){
  downloadButton("saveFig", "Save figure")
}






## Heatmap visualization support options -----
### Heatmap visualization name radio ------
heatmapViz_nameFormat_Radio <- function(){
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
}

### Heatmap TFA color palette -----
heatmapViz_TFAPalette_Select <- function(palettes_edges){
  selectInput(inputId = 'tfa_palette_heatmap', choices = palettes_edges$pal, multiple = FALSE, selected = default_tfa_palette_heatamp,
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
                  )))
}


### Heatmap TFA range -------
heatmapViz_TFARange_Slider <- function(){
  sliderInput(inputId = 'tfa_range_heatmap', value = default_tfa_range, min = default_tfa_min, max = default_tfa_max,
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
                )))
}


### Heatmap expression palette select ------
heatmapViz_expPalette_Select <- function(palettes_edges){
  selectInput(inputId = 'expression_palette_heatmap', choices = palettes_edges$pal, multiple = FALSE, selected = default_expression_heatmap,
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
                )))
}

### Heatmap expression range slider -------
heatmapViz_expRange_Slider <- function(){
  sliderInput(inputId = 'expression_range_heatmap', value = default_expression_range, min = default_expression_min, max = default_expression_max,
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
}


### Heatmap direction radio ------
heatmapViz_direction_Radio <- function(){
  radioButtons(inputId = "heatmap_direction", 
               choices = list("forward" = 1, "backward"= -1), label = tags$div(
                 style = "display: flex; align-items: center;",
                 tags$h4("Heatmap direction"),
                 tags$div(
                   style = "margin-left: 1px;", 
                   bsButton("heatmap_direction_info", "", icon = icon("question-circle", class = "fa-lg"), style = "link"), 
                   bsPopover("heatmap_direction_info", "Additional Info",
                             "Direction of heatmap. Flips colormap high-low axis assignment.",
                             placement = "right",
                             options = list(
                               container = "body",
                               html = TRUE,
                               template = '<div class="popover" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content" style="width: 400px; height: 300px;"></div></div>'
                             )
                   )
                 )))
}


## Heatmap graph visualization options ----
### Heatmap graph viz edge color radio -----
heatmapGraphViz_edgeColor_Radio <- function(){
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
               )))
}

### Heatmap graph viz color range reg ----
heatmapGraphViz_edgeRangeReg_Slider <- function(){
  sliderInput("edge_color_range_heatmap_reg", label = tags$div(
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
      ))), min = -10, max = 10, value = c(-5,5), step = 0.1)
}

### Heatmap Graph viz color range corr ----
heatmapGraphViz_edgeRangeCorr_Slider <- function(){
  sliderInput("edge_color_range_heatmap_corr", label = tags$div(
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
      ))), min = -1, max = 1, value = c(-1,1), step = 0.1)
}


### heatmap Graph edge palette -----
heatmapGraphViz_edgePalette_Select <- function(palettes_edges){
  selectInput(inputId = 'edge_color_palette_heatmap', choices = palettes_edges$pal, multiple = FALSE, selected = default_edge_color_pallette, 
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
  
}




## Heatmap save options ----- 
### Heatmap font size slider -----
heatmapViz_fontSize_Slider <- function(){
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
      ))), min = 0, max = 25, value = 18, step = 1)
}

### Heatmap image height slider -----
heatmapViz_imageHeight_Slider <- function(){
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
                )))
}

### heatmap image width  slider -----
heatmapViz_imageWidth_Slider <- function(){
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
                )))
}

### heatmap file name text ------
heatmapViz_fileName_Text <- function(){
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
                )))
  
}

### heatmap file type select ----
heatmapViz_fileType_Select <- function(){
  selectInput("print_file_type", "File format:",
              choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
              selected = "png")
}

### heaetmap download button ----
heatmapViz_Download <- function()
{
  downloadButton("saveFig2", "Save figure")
}

