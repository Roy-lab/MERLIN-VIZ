# MerlinViz

## Dependencies 
Install these into an R enviroment first!
* shiny
* tools
* tidyverse
* networkD3
* Matrix
* scales
* DT
* webshot
* htmlwidgets

## GUI interface

Users can search the MERLIN network using various features, depending on the selected feature, a meaningful sub network is generated. Feature search include: 

* By Module: A representation of the module is generated with the additional regulated genes.

* By Gene: All neighbors of a gene are shown in addition to any module related to the gene.  (May need to allow searching for regulators of module as often these are not in the module themselves). 

* By GeneList: An extention of the by gene method. A list of genes is uploaded to the server. A subgraph is generated containing all related modules, and neighbors.  (see note above). 

* By GO Term: All modules enriched for a particular GO term are shown. 

* Node Diffusion: A list of genes is uploaded to the server (optionally the file can contain node values such as log expression or q-value). Node diffusion is applied. In the case where no information is provided an arbitrary value of 100 is assigned to each node. The kernel diffusion is applied with selected lambda and and a subgraph with nodes greater than the x% is displayed. x can be set by user. 

In addition to search and node diffusion, the GUI can produce module files similar to the searchForModule function that I have written previously. Specifically when a GeneList is searched, module enrichement is computed for the gene list is computed. 

Finally, an approximate Steiner tree search function has been provided to try and connect unconnect subgraph components when searching for multiple genes. If GO Term or Gene List is used to search network then the user can select to turn on a stiener tree visualization. The Steiner tree is computed and displayed over the subgraph and is displayed with red edges. 

During any visualization, the user can gain additional information about the module and the node by clicking on the a node of interest. This will automatically print the Node Info and Module info to screen. The user can then download this information to file. When more than 1 module is displayed, all module information can be downloaded at once. This would be equivelent to the the searchForModule function that I have been using in the past. 


---

## Notes on implementation and in lab usage: 
Steps to create a GUI:
1. Install library dependencies. 
2. Source aux_functions.R file in R. 
3. run makePostProcessDataStruct (see description below). A list of example files is commented out of the beginning of aux_files. These were used to generate the aspergillus network visualization. This function will create a net_data.Rdata file that will contain the structures required to run the GUI. 
4. run the shiny app with app.R 

The GUI was intended for use by biologist who don't want to interact with a command line interface. However, since all functions used in the GUI are implemented in an independent R file, this file can be sourced into an R command line interface and used to perform all tasks indepedently of the GUI. A list of functions is below:

## Initialization Functions
* **makePostProcessDataStruct(all_nodes_files, edge_list_files, module2gene_files, go_files, module_file, regulator_enrich_file, go_enrich_file, gene2genename_file, gene_desc_file, regulator_list_file, rna_seq_file, atac_seq_file)**: Used to generate two tidyverse data structures that contain all information about the network. The first is a tidygraph containing all merlin edges and nodes, the number of neighbors of each node, a list of neigbhors for each node, the module assignment of node, and the edge confidence. This allows for easy search and use fo the tidygraph algorithm sweet to manipulate data. The second is a Module structure that contains all module information. This is used for efficient lookup of module features.  The function automatically saves this into an Rdata file for later loading. 

  * all_node_files: A list of genes used in inference. I have included genes that do not have any edges as they may be in a module. All genes in the orginal structure inference step. 

  * edge_list_files: A list of merlin edges. Format: <source gene from all_node_file> <target gene from all_node_file> <weight> using tab delimiters. 

  * module2gene_file: a node by node assignment to modules: Module assignment file from merlin. Format <all_node_file gene> <module assigment> using tab delimiter. 

  * go_file: A list of all GO terms assocatiated with genes. Format <Gene name matching all_nodes_file> <GOTerm> <TermLevel> using a tab delimiter. There is a header on this file. 

  * module_file: A file containing a list of genes per module. Should have the module number followed by a list of genes. Format: Cluster<module id> <gene lists using # delimiter>. 

  * regulator_enrich_file: The regulator enrichment of MERLIN modules. Format from enrich analyzer using the module_file.

  * go_enrich_file: The go enrichment of MERLIN modules. Format from enrich anaylzer using the module_file. 

  * gene2genename_file: This file contains two columns. The first column is contains a formal gene name. The second file contains a common gene name. This has use in some plant and fungal species. In Aspergillus, all names have a formal id AFUA_#G#####, which is useful for search but less recognizable then a common name. If you would like to ignore this feature just include a file with the same names in two columns (human etc). 

  * gene_desc_file: This file contains a one-line discription of each gene. e.g. RING finger protein. Format should be <formal gene name from all nodes file> \tab Description (Cannnot include tabs)

  * regulator_list_file: A file containing a list of genes from all_node_file that are designated as regulators in the model specification. Format is a single column of genes.  
  
  * rna_seq_file: A file containing a list of feature points and gene names. The first column should be labeled Names, and have features that match the genes in all_nodes file. The remaining columns must be labeled and presented in the order of the desired display.
 
 * atac_seq_file: A file containing a list of features points and gene_names. The first column should be labeled Names, and have features that match teh genes in all_nodes. The remaining columns must be labeled and presented in the order of the desired display. 
 
* **makeLaplacian(Net)**: Generates the laplacian of the graph. Requires the Net object from makePostProcessDataStruct.

* **MakeKernel(L, lambda)**: Generates a diffusion kernel from laplacian with hyperparameter lambda. These can be saved and loaded for later use.

### Search functions

* **searchForModules(Module, moduleID)**: Returns a list of genes that are related to a module. Thes include the regulators of the module and genes in the module.

* **searchForGenes(Net, Module, gene)**: Returns a list of genes related to a gene of interest. This inlcudes genes within the neighborhood of the gene, all module genes. 

* **searchForGeneList(Net, Module, gene_list)**: Return the search for genes for a list of genes given by gene_list. 

* **computeEnrichment(Module, gl, num_genes)**: Computes the module enrichment for a list of genes gl. The number of genes in gl also needs to be provided. The enrichment is computed via the hypergeometric test and corrected via bonferonni correction. 

# Diffusion Functions

* **generateScoreVector(Net, gene_list)**: Takes a gene list and generates an automatic score vector used in diffusion. All gene are assigned the arbitrary value of 100 if they are on the gene list. All other genes are assigned 0.

* **loadScoreVector(Net, score_data)**: Similar in usage to the generateScoreVector, loadScoreVector allows users to submit some scores of interested included q-value or abs(log(fold_change)). Returns a vector with the scores given by the input file, otherwise assigns 0. 

* **computeDiffusionScore(Net, score_data, kernel)**: The diffusion scores are added to the Net structure after node diffusion using the given input kernel.]

###Wrappers for subgraph generation using search for methods.

Each of these produces a tidygraph subgraph that is generated using the corresponding search for function:

* **moduleSubgraph(Net, Module, module_id)**:  Makes a subgraph of the specific module given by module id. 

* **geneSubgraph(net, Module, gene_list)**: Makes a subgraph using the search for genes function.

* **geneListSubgraph(net, Module, gene_list)**: Makes a subgraph using the search for genes_list function.

* **goSubgraph(Net, Module, enrich_2_module, go_term)**: Makes a subgraph using the the given go terms. 

* **diffScoreSubgraph(Net, Percentile)**: Uses diffusion scores with cutoff percentile to generate a subgraph. 

* **induceSubgraph(Net, list)**: A helper function that induces a subgraph given a list of genes. 

* **graph2NodeEdgeTable(Net)**: Takes in a tidygraph network and produces tibbles (tidyverse table structure) containing node and edge lists. Used to generate required format for graphing. Note that the edge indexing is change for 1 to 0 indexing for usage with NetworkD3 display module. 

### Steiner Tree functions

* **getDistMatrix(Net, gene_list)**: given a gene list, gives the distance form each node in the list to every other node in the network. Uses a breadth first search method. This is used to generate psuedo_stiener trees. 

* **buildStienerTrees(Net, gene_list)**: A wrapper function for generation of the Steiner tree. The dist matrix is generated, then the closest 2 nodes are connected. If more than 1 pair has the smallest distance the first in numerical order by index is selected. Nodes are then successively added to the tree by connected the next closest node. This continues until all nodes are as connected as possible. Note this function no longer requires only the largest connected component. Instead, a forest is generated if there are two unconnected components of interest. Assumes undirected graph structure. 

### PrinterFunctions

* **printNodeInfo(Net, node_name)**: Prints information of Node in HTML format. Information includes name, module, associated GO term, and neighbors. 

* **printModuleInfo(Module, module_id, gene_list, genes)**: Prints information about module in HTML format. Information includes module id, module enrichment if computed, module genes, enriched regulators, and module enriched go term. 

* **PrintAllModuleInfo(Subnet, Module, gene_list, genes)**: Prints a list of all modules in HTML format. This is used when a subgraph contains more than one module and makes successive callse to printModuleInfo. 

* **GetModuleID(Net, node_name)**: returns the module id of the gene given by node name. 



