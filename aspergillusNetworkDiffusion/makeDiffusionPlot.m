function [p, S]=makeDiffusionPlot(head_dir, scores_file, common_names_file,  cut_off_q, topk, outprefix)
node_index=readtable(sprintf('%s/aspergillus_network_I02_node_index.txt', head_dir), 'ReadVariableNames', 0);
edge_list=importdata(sprintf('%s/aspergillus_network_I02_edge_list.txt', head_dir));
common_names=readtable(common_names_file, 'Delimiter', 'tab', 'ReadVariableNames', 0);

common_names_map=containers.Map(common_names.Var1, common_names.Var2);

node_index=table2cell(node_index);
for i=1:length(node_index)
    if common_names_map.isKey(node_index{i,2});
        node_index(i,2)={common_names_map(node_index{i,2})};
    end
end

node_index_map=containers.Map(node_index(:,1), node_index(:,2));


G=graph(edge_list(:,1), edge_list(:,2));

diffusion_score=readtable(scores_file, 'Delimiter', 'tab', 'ReadVariableNames', 1);
diffusion_score=table2array(diffusion_score);

Name={};
for i=1:height(G.Nodes)
    Name(i)={strrep(node_index_map(i), '_', '\_')};
    idx=find(diffusion_score(:,1)==i);
    G.Nodes.Score(i)=diffusion_score(idx,3);
end
G.Nodes.Name=Name';
cut_off=quantile(G.Nodes.Score, cut_off_q);
Gidx=find(G.Nodes.Score>=cut_off);
S=subgraph(G,Gidx);
%p=plot(S, 'Layout', "force", 'MarkerSize', 5);
p=plot(S, 'MarkerSize', 3);
layout(p, 'force', 'UseGravity', true);
p.NodeCData=S.Nodes.Score;
[~, top]=sort(S.Nodes.Score, 'descend');
saveas(gcf, sprintf('%s.fig', outprefix));
labelnode(p, top(1:topk), S.Nodes.Name(top(1:topk)));
saveas(gcf, sprintf('%s.png', outprefix));
