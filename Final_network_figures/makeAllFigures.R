#### Make Figures for Asp Paper
library(Matrix)
## Clear everything
rm(list = ls())
## Load in network and functions 
source('aux_functions.R')
source('printerFunction.R')
load('net_data.Rdata')

firstUpper <- function(x){
  substring(x, 1, 1) <- toupper(substring(x,1,1))
}


###Figure 3.b
subNet <- geneListSubgraph(Net, Module, c('AFUA_2G01260'), search_additional = c('mod'))
print_disp_names <- c('cyp51A', 'erG25B', 'hyd1', 'srbA', 'srbB', 'erG3', 'erG25', 'fhpA', 'erG1', 
                      'hem13', 'niiA', 'AFUA_5G06120_nca', 'AFUA_3G12190', 'bna4', 'srb5', 'hem14', 'exG4', 
                      'erG3A', 'pre4', 'AFUA_7G04740', 'AFUA_6G02180')


subNet <- subNet %N>% mutate(component = group_components()) 
keep_component <- subNet %N>% as_tibble() %>% 
  group_by(component) %>% 
  summarise(count  = n()) %>% 
  filter(count >= 2) %>% 
  pull(component)
subNet <- subNet %N>% filter(component %in% keep_component )
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))

gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'geneSuper', 
                edge_color_by = NA, node_color_palette = 'Dark2', 
                node_size_by = NA, max_node_size = 4, 
                layout = 'dh', focus_nodes = list(), 
                font_size = 4, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                expand_x = 1.5, expand_y = 0)

ggsave('fig3_srbA_srbB_white_font.png', gg, width = 8.5, height = 5, dpi = 900)

gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'geneSuper', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = NA, max_node_size = 4, 
                      layout = 'dh', focus_nodes = list(), 
                      font_size = 6, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0, font_color = '#000000', unlab_color = '#dddddd')


ggsave('fig3_srbA_srbB_black_font.png', gg, width = 8.5, height = 5, dpi = 900)


### Figure 4.A 
subNet <- geneListSubgraph(Net, Module, c('AFUA_1G10080'), search_additional = c('neigh'))
print_disp_names <- c('zrfA', 'zrfC', 'zafA', 'sod3', 'zrfB', 'AFUA_5G02010', 'AFUA_1G14700')

subNet <- subNet %N>% mutate(component = group_components()) 
keep_component <- subNet %N>% as_tibble() %>% 
  group_by(component) %>% 
  summarise(count  = n()) %>% 
  filter(count >= 2) %>% 
  pull(component)
subNet <- subNet %N>% filter(component %in% keep_component )
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))

gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'geneSuper', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = NA, max_node_size = 4, 
                      layout = 'dh', focus_nodes = list(), 
                      font_size = 6, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig4a_zrfA_white_font.png', gg, width = 7, height = 4, dpi = 900)


### Fig 4.B
subNet <- geneListSubgraph(Net, Module, c('AFUA_5G03920'), search_additional = c('mod'))
print_disp_names <- c('sidA', 'sidD', 'sidH', 'sidF', 'sidJ', 'sidG', 'mdr4', 'mirB', 'amcA', 'sidI', 'sit1', 'fre2', 'sit2', 
                      'mirD', 'hsf1_nca', 'cccA', 'ccp1', 'estB', 'hapX', 'sitT')

subNet <- subNet %N>% mutate(component = group_components()) 
keep_component <- subNet %N>% as_tibble() %>% 
  group_by(component) %>% 
  summarise(count  = n()) %>% 
  filter(count >= 2) %>% 
  pull(component)
subNet <- subNet %N>% filter(component %in% keep_component )
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))
gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'geneSuper', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = NA, max_node_size = 4, 
                      layout = 'dh', focus_nodes = list(), 
                      font_size = 6, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig4b_hapX_white_font.png', gg, width = 8, height = 5, dpi = 900)


### Fig 4.c
subNet <- geneListSubgraph(Net, Module, c('AFUA_8G00420'), search_additional = c('mod')) 
print_disp_names <- c('fmaB', 'Fma-KR', 'fsqG', 'fsqF', 'psoB', 'psoF', 'hsf1_nca', 'fmaG', 'fmaD', 'fmaC', 'fmaF', 'fapR', 
                      'psoE', 'psoC', 'fmaA', 'fsqC', 'fsqB', 'psoD', 'nrps14', 'nsdD_nca', 'lreA_nca')
subNet <- subNet %N>% mutate(component = group_components()) 
keep_component <- subNet %N>% as_tibble() %>% 
  group_by(component) %>% 
  summarise(count  = n()) %>% 
  filter(count >= 2) %>% 
  pull(component)
subNet <- subNet %N>% filter(component %in% keep_component )
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' ')) %>%
  mutate(display_name = ifelse(display_name == 'Fma-KR', 'fma-KR', display_name))
gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'geneSuper', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = NA, max_node_size = 4, 
                      layout = 'kk', focus_nodes = list(), 
                      font_size = 3, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig4c_SM_white_font.png', gg, width = 10.5, height = 4, dpi = 900)

### Fig 5b 
subNet <- geneListSubgraph(Net, Module, c('AFUA_3G11990_nca'), search_additional = c('mod')) 
print_disp_names <- c('GliZ', 'nrps9', 'GliI', 'GliH', 'GliP', 'GliC', 'GliM', 'GliG', 'GliK', 'GliN', 'GliK', 'GliN',
                      'GliA', 'GliT', 'GliF', 'GtmA', 'hsf1_nca', 'AFUA_3G11990_nca')
subNet <- subNet %N>% mutate(component = group_components()) 
keep_component <- subNet %N>% as_tibble() %>% 
  group_by(component) %>% 
  summarise(count  = n()) %>% 
  filter(count >= 2) %>% 
  pull(component)
subNet <- subNet %N>% filter(component %in% keep_component )
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))
gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'geneSuper', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = NA, max_node_size = 4, 
                      layout = 'kk', focus_nodes = list(), 
                      font_size = 6, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig5b_gli_white_font.png', gg, width = 6, height = 6, dpi = 900)


### Fig 6a 
score_list <- read_csv(file = 'Afumigatus_sLCO_DEGs_30minutes_edgeR_read_in_logFC.csv', col_names = c("feature", "score"))
load('k10.Rdata')
Net <- computeDiffusionScore(Net, score_list, k10_sparse)
min_targets <- 5 
top_regs <- 10
subNet <- diffScoreSubgraph(Net, min_targets, top_regs)
nodes <-  Net %N>% as_tibble() %>%  arrange(desc(score)) %>% filter(regulator == 'scr') %>% filter(degree >= min_targets)
print_disp_names <- nodes  %>% slice(1:top_regs) %>% pull(`Common Name`)
subNet <- subNet %N>% mutate(module = as.character(module)) %>% mutate(module = str_replace(module, '-9999', 'Unlabeled'))
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))
gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'module', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = 'score', max_node_size = 10, 
                      layout = 'dh', focus_nodes = list(), 
                      font_size = 8, nudge_y = 0, text_angle  = 0, show_legend = TRUE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig6a_lco_30min_abs_log_FC_white_font_legend.png', gg, width = 15, height = 10, dpi = 900)

gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'module', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = 'score', max_node_size = 10, 
                      layout = 'dh', focus_nodes = list(), 
                      font_size = 8, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig6a_lco_30min_abs_log_FC_white_font.png', gg, width = 15, height = 10, dpi = 900)



### Fig 6B 
score_list <- read_csv(file = 'Afumigatus_sLCO_DEGs_2hours_edgeR _read_in_logFC.csv', col_names = c("feature", "score"))
load('k10.Rdata')
Net <- computeDiffusionScore(Net, score_list, k10_sparse)
min_targets <- 5 
top_regs <- 10
subNet <- diffScoreSubgraph(Net, min_targets, top_regs)
nodes <-  Net %N>% as_tibble() %>%  arrange(desc(score)) %>% filter(regulator == 'scr') %>% filter(degree >= min_targets)
print_disp_names <- c( unlist(nodes  %>% slice(1:top_regs) %>% pull(`Common Name`)), 'atfA_nca')
subNet <- subNet %N>% mutate(module = as.character(module)) %>% mutate(module = str_replace(module, '-9999', 'Unlabeled'))
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))
gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'module', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = 'score', max_node_size = 10, 
                      layout = 'dh', focus_nodes = list(), 
                      font_size = 8, nudge_y = 0, text_angle  = 0, show_legend = TRUE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig6b_lco_2hrs_abs_log_FC_white_font_legend.png', gg, width = 11, height = 10, dpi = 900)

gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'module', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = 'score', max_node_size = 10, 
                      layout = 'dh', focus_nodes = list(), 
                      font_size = 8, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig6b_lco_2hrs_abs_log_FC_white_font.png', gg, width = 10, height = 10, dpi = 900)





### Fig 7a 

subNet <- geneListSubgraph(Net, Module, c('AFUA_3G14420','AFUA_2G01870'), search_additional = c('mod')) 
print_disp_names <- c('chsA', 'chsG', 'rfeF', 'Gfa1', 'rGd1', 'GpaA', 'medA', 'AFUA_7G04300', 'AFUA_4G13600', 'gpaA')
subNet <- subNet %N>% mutate(component = group_components()) 
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))
subNet <- subNet %N>% mutate(module = as.character(module)) %>% mutate(module = str_replace(module, '-9999', 'Unlabeled'))
gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'module', 
                      edge_color_by = NA, node_color_palette = 'Dark2', 
                      node_size_by = NA, max_node_size = 4, 
                      layout = 'kk', focus_nodes = list(), 
                      font_size = 6, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)

ggsave('fig7a_no_steiner_white_font.png', gg, width = 6, height = 6, dpi = 900)


### Fig 7b
subNet <- geneListSubgraph(Net, Module, c('AFUA_3G14420','AFUA_2G01870'), search_additional = c('mod', 'stein')) 
print_disp_names <- c('chsA', 'chsG', 'rfeF', 'gfa1', 'rGd1', 'gpaA', 'medA', 'AFUA_7G04300', 'AFUA_4G13600', 'gpaA')
subNet <- subNet %N>% mutate(component = group_components()) 
keep_component <- subNet %N>% as_tibble() %>% 
  group_by(component) %>% 
  summarise(count  = n()) %>% 
  filter(count >= 2) %>% 
  pull(component)
subNet <- subNet %N>% filter(component %in% keep_component )
subNet <-subNet %N>% mutate(display_name = ifelse(`Common Name` %in% print_disp_names, `Common Name`, NA)) %>% mutate(display_name = str_replace(display_name, '_', ' '))
subNet <- subNet %N>% mutate(module = as.character(module)) %>% mutate(module = str_replace(module, '-9999', 'Unlabeled'))
gg <- makeSubNetGraph(subNet, names_in_nodes = TRUE, node_color_by = 'module', 
                      edge_color_by = 'is_steiner', node_color_palette = 'Dark2', 
                      node_size_by = NA, max_node_size = 4, 
                      layout ='nicely', focus_nodes = list(), 
                      font_size = 6, nudge_y = 0, text_angle  = 0, show_legend = FALSE, 
                      expand_x = 1.5, expand_y = 0)
ggsave('fig7b_steiner_white_font.png', gg, width = 6, height = 6, dpi = 900)


