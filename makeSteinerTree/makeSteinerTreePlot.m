function fig=makeSteinerTreePlot(infile, name_map, outprefix)
T=readtable(infile, 'Delimiter', 'tab', 'ReadVariableNames', 0);
T=table2cell(T);

common_names=readtable(name_map, 'Delimiter', 'tab', 'ReadVariableNames', 0);
name_map=containers.Map(common_names.Var1, common_names.Var2);


for i=1:length(T)
    if contains(T{i,1}, '_nca')
        temp=strrep(T{i,1}, '_nca', '');
        if name_map.isKey(temp)
           T(i,1)={name_map(temp)};
        end
        T(i,1)={strcat(temp, '_nca')};
    else
        if name_map.isKey(T{i, 1})
           T(i,1)={name_map(T{i,1})};
        end
    end
    
    if contains(T{i,2}, '_nca')
        temp=strrep(T{i,1}, '_nca', '');
        if name_map.isKey(temp)
           T(i,2)={name_map(temp)};
        end
        T(i,2)={strcat(temp, '_nca')};
    else
        if name_map.isKey(T{i, 2})
           T(i,2)={name_map(T{i,2})};
        end
    end
    T(i,1)={strrep(T{i,1}, '_', '\_')};
    T(i,2)={strrep(T{i,2}, '_', '\_')};    
end


G=graph(T(:,1), T(:,2));
fig=plot(G)
layout(fig, 'layered', 'Direction', 'down');
saveas(gcf, sprintf('%s.fig', outprefix));

labelnode(fig,G.Nodes.Name, G.Nodes.Name)
fig.NodeFontSize=8;
ax=gca;
ax.Position=[0 0 1 1];
saveas(gcf, sprintf('%s.png', outprefix));