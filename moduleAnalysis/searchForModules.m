function [Modules]=searchForModules(Gene_file, Name_Mapping, common_name_in, common_name_out)
%%% This function is used to search a list of modules from Merlin for
%%% specific genes. Four files are required to be in the directory of this
%%% function. 
%%%
%%% Required file:
%%% consensus_module_0_3_geneset_enrichAnalyzer.xlsx - Contains merlin
%%% clusters. The first column contains the cluster id. The second column
%%% contains the genes within the module. Genes are deliminated by ', '
%%% 
%%% merlin.0_8.0_3.details.xlsx - This file contains the genes that are
%%% enriched regulators of the cluster. The first column contains the
%%% module for interest the second column contains the Gene ID of the
%%% Regulator. The last column contains the target genes of the enriched
%%% regulator. 
%%%
%%% go.0_3.details.xlsx - This files contains the enriched go term of
%%% modules 
%%% 
%%% all_nodes.txt - A list of all TF, possible targets, and TFAs. 
%%%
%%% Input:
%%% Gene_file - The path to a file containing a list of genes of interest.
%%% These could be DE genes or specific genes. The file should be delimited
%%% by tab or be in csv format. Only the first column is considered and
%%% gene names need to be in the following format AFUA_#G#####
%%%
%%% Outputs:
%%% Modules - A matlab structure that contains the data storted in the
%%% three different required files organized by structure. If no enrichment
%%% is found for module the subsiquential fields are empty.
%%%
%%% <Gene_file>_module_summary.txt - This is a formatted file containing
%%% information about the genes in gene_file. The first row contains all 
%%% modules that contain genes from gene_file sorted by the p-value of the 
%%% fdr corrected hypergeometic test. Then for each module containing a
%%% gene on the gene list summary information is printed, including the
%%% enrichment and GO term. 

%% Read in file info:
head_dir='/mnt/dv/wid/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100'
%head_dir='/Volumes/wid-1/projects7/Roy-Aspergillus/Results/RnaSeq/MerlinpResults/Afum_RnaSeq_results/Condor_results/PostBatchCorr_I02/Lambda_0100';
Clusts=readtable(sprintf('%s/consensus_module_0_3_geneset_enrichAnalyzer.txt', head_dir), 'ReadVariableNames', 0);
Clusts.Properties.VariableNames={'ModuleID', 'GeneInTheModule'};

Merlin_info=readtable(sprintf('%s/Enrichments/merlin.0_8.0_3_details.txt', head_dir), 'ReadVariableNames', 0);
T=table2cell(Merlin_info(:,10));
for i=1:length(T)
    T(i)={split(T(i), ';')};
end
Merlin_info(:,10)=T;
Merlin_info.Properties.VariableNames(1)={'ModuleID'};

GO_info=readtable(sprintf('%s/Enrichments/go.0_3_details.txt', head_dir));
T=table2cell(GO_info(:,10));
for i=1:length(T)
    T(i)={split(T(i), ';')};
end
GO_info(:,10)=T;
GO_info.Properties.VariableNames(1)={'ModuleID'};
GO_info.Properties.VariableNames(2)={'GOTerm'};



Common_name=readtable(Name_Mapping, 'Delimiter', 'tab', 'ReadVariableNames', 0);
comm_2_gen=containers.Map(Common_name.Var2, Common_name.Var1);
gen_2_comm=containers.Map(Common_name.Var1, Common_name.Var2);

Links2Labeled=containers.Map('KeyType', 'char', 'ValueType', 'char');
Links2Standard=containers.Map('KeyType', 'char', 'ValueType', 'char');
Links_list=table2cell(readtable('I02_required_files/listing_of_folder_Hmap_out.0_8.0_3.xlsx'));
for i=1:length(Links_list)
    if contains(Links_list(i, 1), 'labels.svg')
        s=split(Links_list(i,1), '_');
        clust=s{1};
        Links2Labeled(clust)=Links_list{i,2};
    elseif contains(Links_list(i,1), '.svg')
        s=split(Links_list(i,1), '.');
        clust=s{1};
        Links2Standard(clust)=Links_list{i,2}; 
    end
end





Gene_list_data=readtable(Gene_file, 'Delimiter', 'tab', 'ReadVariableNames', 0);
Gene_list=table2cell(Gene_list_data(:,1));

if common_name_in
    for i=1:length(Gene_list)
        if isKey(comm_2_gen, Gene_list{i})
            Gene_list{i}=comm_2_gen(Gene_list{i});
        else
            sprintf('Warning: general name for %s is not in map file', Gene_list{i})
        end
        
    end  
end

all_nodes=importdata(sprintf('%s/all_nodes_0_8.txt',head_dir));
nca_nodes=contains(all_nodes,'nca');
nca_split=split(all_nodes(nca_nodes),'_');
nca_genes=strcat(nca_split(:,1),'_',nca_split(:,2));
add_nca=intersect(Gene_list,nca_genes);
Gene_list(end+1:end+length(add_nca))=strcat(add_nca,'_nca');

Gene_list=intersect(Gene_list, all_nodes);



%% Generate Module structure
Modules=struct([]);
for i=1:height(Clusts)
    Modules(i).ID=(Clusts.ModuleID(i));  % Module id 
    Modules(i).Genes=split(Clusts.GeneInTheModule(i),'#'); % cell array of Gene name strings
    
    Merlinidx=find(contains(Merlin_info.ModuleID, Modules(i).ID)); % find all enriched regulators
    if ~isempty(Merlinidx)
    Modules(i).EnrichedRegulators=table2cell(Merlin_info(Merlinidx,2)); % cell array of enriched regulators
    Modules(i).regulator_pvalue=table2array(Merlin_info(Merlinidx,3));     % matrix of associated pvalues, ordered by regulators
    Modules(i).regulator_corrected_pvalue=table2array(Merlin_info(Merlinidx,4)); %Corrected p-values
    Modules(i).regulator_enrichment=table2array(Merlin_info(Merlinidx,9)); % Fold Enrichment
    Modules(i).regulated_genes=table2cell(Merlin_info(Merlinidx,10));      % regulated genes
    end
    
    GOidx=find(contains(GO_info.ModuleID, Modules(i).ID));
    if ~isempty(GOidx)
        Modules(i).GOTerms=GO_info.GOTerm(GOidx);           %List of GO terms
        Modules(i).GO_pvalue=table2array(GO_info(GOidx,3)); % Associated pvalue
        Modules(i).GO_corrected_pvalue=table2array(GO_info(GOidx,4)); % corrected p value
        Modules(i).GO_Enrichment=table2array(GO_info(GOidx,9)); % Fild enrichment
        Modules(i).GO_Genes=(table2cell(GO_info(GOidx,10)));  % Genes corresponding to Go.
    end
end


%% Find Modules that contain genes from genelist. Clustidx is a 1-hot matrix for genes in modules. 
Clustidx=zeros(length(Gene_list),length(Modules));
for i=1:length(Gene_list)
    for j=1:length(Modules)
        if ~isempty(Modules(j).EnrichedRegulators)
            Clustidx(i,j)=sum(contains(Modules(j).Genes, Gene_list{i}))+sum(contains(Modules(j).EnrichedRegulators, Gene_list{i}));
        else
            Clustidx(i,j)=sum(contains(Modules(j).Genes, Gene_list{i}));
        end
    end
end


[idxI,idxJ]=find(Clustidx);
unique_clusts=unique(idxJ);          %List of unique modules
%% Compute Hypergeometric test value
tot_nodes=length(all_nodes);
num_search_genes=length(Gene_list);

num_genes_list=zeros(1,length(unique_clusts));
gene_list_enrichment_pvalue=zeros(1,length(unique_clusts));
gene_list_fold_change=zeros(1, length(unique_clusts));
clust_size=zeros(1, length(unique_clusts));
for i=1:length(unique_clusts)
    num_genes_list(i)=sum(unique_clusts(i)==idxJ); %Number of gene from gene list in modules
    clust_size(i)=length(Modules(unique_clusts(i)).Genes);
    gene_list_enrichment_pvalue(i)=hygecdf(num_genes_list(i)-1, tot_nodes, num_search_genes, clust_size(i),'upper');
    gene_list_fold_change(i)=(num_genes_list(i)/clust_size(i))/(num_search_genes/tot_nodes);
end
[gene_list_enrichment_corrected_pvalue]=mafdr(gene_list_enrichment_pvalue,'BHFDR',1);
[~, sort_idx]=sort(gene_list_enrichment_corrected_pvalue,'ascend'); 

for i=1:length(unique_clusts)
    Modules(unique_clusts(i)).Gene_List_Enrichment=gene_list_fold_change(i);
    Modules(unique_clusts(i)).gene_list_pvalue=gene_list_enrichment_pvalue(i);
    Modules(unique_clusts(i)).gene_list_corrected_pvalue=gene_list_enrichment_corrected_pvalue(i);
end


%% Change Gene Names
if common_name_out
    for i=1:length(Modules)
        %%%% Change Gene Names %%%%
        for j=1:length(Modules(i).Genes)
            gene_name=Modules(i).Genes{j};
            if contains(gene_name, 'nca')
                gs=split(gene_name,'_');
                gene_name=strjoin(gs(1:2), '_');
                if isKey(gen_2_comm, gene_name)
                    Modules(i).Genes{j}=strcat(gen_2_comm(gene_name), '_nca');
                end
            else
                if isKey(gen_2_comm, gene_name)
                    Modules(i).Genes{j}=gen_2_comm(gene_name);
                end
            end
        end
        %%% Change Regulator Names %%%%%
        for j=1:length(Modules(i).EnrichedRegulators)
            gene_name=Modules(i).EnrichedRegulators{j};
            if contains(gene_name, 'nca')
                gs=split(gene_name,'_');
                gene_name=strjoin(gs(1:2), '_');
                if isKey(gen_2_comm, gene_name)
                    Modules(i).EnrichedRegulators{j}=strcat(gen_2_comm(gene_name), '_nca');
                end
            else
                if isKey(gen_2_comm, gene_name)
                    Modules(i).EnrichedRegulators{j}=gen_2_comm(gene_name);
                end
            end
        end
       %%%% Change Regulated Names
       for j=1:length(Modules(i).regulated_genes)
           Gene_set=Modules(i).regulated_genes{j};
           for k=1:length(Gene_set)
               gene_name=Gene_set{k};
               if contains(gene_name, 'nca')
                    gs=split(gene_name,'_');
                    gene_name=strjoin(gs(1:2), '_');
                    if isKey(gen_2_comm, gene_name)
                        Modules(i).regulated_genes{j}{k}=strcat(gen_2_comm(gene_name), '_nca');
                    end
               else
                    if isKey(gen_2_comm, gene_name)
                        Modules(i).regulated_genes{j}{k}=gen_2_comm(gene_name);
                    end
               end
           end
       end
       %%% Change GO Names
       for j=1:length(Modules(i).GO_Genes)
           Gene_set=Modules(i).GO_Genes{j};
           for k=1:length(Gene_set)
               gene_name=Gene_set{k};
               if contains(gene_name, 'nca')
                    gs=split(gene_name,'_');
                    gene_name=strjoin(gs(1:2), '_');
                    if isKey(gen_2_comm, gene_name)
                        Modules(i).GO_Genes{j}{k}=strcat(gen_2_comm(gene_name), '_nca');
                    end
               else
                    if isKey(gen_2_comm, gene_name)
                        Modules(i).GO_Genes{j}{k}=gen_2_comm(gene_name);
                    end
               end
           end
       end
       
    end
    
    for i=1:length(Gene_list)
        gene_name=Gene_list{i};
        if contains(gene_name, 'nca')
            gs=split(gene_name,'_');
            gene_name=strjoin(gs(1:2), '_');
            if isKey(gen_2_comm, gene_name)
               Gene_list{i}=strcat(gen_2_comm(gene_name), '_nca');
            end
        else
            if isKey(gen_2_comm, gene_name)
                Gene_list{i}=gen_2_comm(gene_name);
            end
        end
    end 
end



%% Write Formatted Summary
file_path=sprintf('%s/moduleAnalysis/', head_dir);
mkdir(file_path)

[~, basename, ~]=fileparts(Gene_file);
if common_name_out
    outfile=sprintf('%s/%s', file_path, strcat(basename,'_common_name_module_summary.txt'));
    %outfile2=sprintf('%s/%s', file_path, strcat(basename, '_common_name_module_genes.txt'));
    outfile3=sprintf('%s/%s', file_path, strcat(basename, '_sig_modules.txt'));
else
    outfile=sprintf('%s/%s', file_path, strcat(basename,'_module_summary.txt'));
    %outfile2=sprintf('%s/%s', file_path, strcat(basename,'_module_genes.txt'));
    outfile3=sprintf('%s/%s', file_path, strcat(basename, '_sig_modules.txt'));
end

delete(outfile);
fileID=fopen(outfile, 'w');
fprintf(fileID, "Module Summary Report \n");

%for i=1:length(num_genes_list)
%   id= Modules(unique_clusts(sort_idx(i))).ID;
%   fprintf(fileID,"%s:%d\t", id{1}, gene_list_enrichment_corrected_pvalue(sort_idx(i)));
%end
tot_genes={};
%fprintf(fileID,'\n\n\n');
%fprintf(fileID,"Begin Module Details \n");

fprintf(fileID,"___________________________________________________________________________________________________\n");
for i=1:length(num_genes_list)
    idx=unique_clusts(sort_idx(i));
    id= strjoin(Modules(idx).ID);
    gene_list_enrich=Modules(idx).Gene_List_Enrichment;
    gene_list_pval=Modules(idx).gene_list_pvalue;
    gene_list_corr_pval=Modules(idx).gene_list_corrected_pvalue;
    matching_genes = strjoin(Gene_list(idxI(idxJ==idx)),', ');
    all_genes = strjoin(Modules(idx).Genes, ', ');
    tot_genes=[tot_genes; Modules(idx).Genes];
    fprintf(fileID, 'Cluster id: %s \nEnrichment with Genes from Gene List:%.02f \t p-value:%d \t Corrected p-value:%d \nGenes from Gene List: %s \nAll Genes: %s \n',...
        id, gene_list_enrich, gene_list_pval, gene_list_corr_pval, matching_genes, all_genes);
    fprintf(fileID, '\nRegulators:');
    if ~isempty(Modules(idx).EnrichedRegulators)
        for j=1:length(Modules(idx).EnrichedRegulators)
            Regulator=Modules(idx).EnrichedRegulators{j};
            if contains(Regulator,Gene_list)
                isOnGeneList='Yes';
            else
                isOnGeneList='No';
            end
            regulator_pvalue=Modules(idx).regulator_pvalue(j);
            regulator_corr_pvalue=Modules(idx).regulator_corrected_pvalue(j);
            regulator_enrichment=Modules(idx).regulator_enrichment(j);
            target_genes_on_list=strjoin(intersect(Modules(idx).regulated_genes{j},Gene_list),', ');
            all_target_genes=strjoin(Modules(idx).regulated_genes{j},', ');
            
            fprintf(fileID, '\nRegulator: %s \t Is On Gene List: %s \t Regulator Enrichment: %.02f \t p-value:%d \t Corrected p-value:%d \n',...
                Regulator,isOnGeneList, regulator_enrichment, regulator_pvalue, regulator_corr_pvalue);
            fprintf(fileID, '\tTarget Genes on Gene List: %s \n', target_genes_on_list);
            fprintf(fileID, '\tAll Target Genes: %s \n\n', all_target_genes);
        end
    else
        fprintf(fileID, ' No Regulators \n');
    end
    
    fprintf(fileID, 'GOTerms: ');
    if ~isempty(Modules(idx).GOTerms)
        for j=1:length(Modules(idx).GOTerms)
            GO_Term=Modules(idx).GOTerms{j};
            GO_pvalue=Modules(idx).GO_pvalue(j);
            GO_corr_pvalue=Modules(idx).GO_corrected_pvalue(j);
            GO_enrich=Modules(idx).GO_Enrichment(j);
            GO_genes_on_list=strjoin(intersect(Modules(idx).GO_Genes{j},Gene_list),', ');
            GO_all_genes=strjoin(Modules(idx).GO_Genes{j},', ');
            fprintf(fileID, '\n GO: %s \t GO Enrichment: %.02f \t p-value:%d \t Corrected p-value:%d \n',...
                GO_Term, GO_enrich, GO_pvalue, GO_corr_pvalue);
            fprintf(fileID, '\tGO Annotated Genes from Gene List: %s \n', GO_genes_on_list);
            fprintf(fileID, '\tAll GO Annotated Genes: %s \n\n', GO_all_genes);
        end
    else
       fprintf(fileID, 'No GO');
    end
    fprintf(fileID, '\n\n');
    if Links2Labeled.isKey(id)
    fprintf(fileID, 'Links to images: \n');
    fprintf(fileID, 'Labeled: %s \n', Links2Labeled(id));
    fprintf(fileID, 'Reference Number: %s \n', Links2Standard(id));
    fprintf(fileID,'\n\n');
    end
    fprintf(fileID,"___________________________________________________________________________________________________\n");
end




fclose(fileID);

sig_modules=[];
for i=1:length(Modules)
if Modules(i).gene_list_corrected_pvalue<.05
sig_modules(end+1)=i;
end
end

sc=[Modules(sig_modules).ID]';
sc=cell2table(sc);
writetable(sc, outfile3, 'WriteVariableNames',0,'Delimiter', 'tab'); 

tot_genes=cell2table(tot_genes);
%writetable(tot_genes,outfile2,'WriteVariableNames',0,'Delimiter', 'tab')
