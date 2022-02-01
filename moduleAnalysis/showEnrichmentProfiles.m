function showEnrichmentProfiles(files, names, avg_module_hmap)
T=importdata(avg_module_hmap);
C=split(T);
clusts=split(C(:,1), '||');
unique_clusts=unique(clusts(:,1));
ID=split(unique_clusts, 'r');
ID=str2num(cell2mat(ID));

row_map=containers.Map(unique_clusts, 1:length(unique_clusts));
sig_Assign=zeros(length(unique_clusts), length(files));

for i=1:length(files)
    sig_clusts=importdata(files{i});
    for j=1:length(sig_clusts)
        if  row_map.isKey(sig_clusts{j})
            idx=row_map(sig_clusts{j});
            sig_Assign(idx, i)=1;
        end
    end
end

of=fopen(outfile,'w');

for i=1:length(unique_clusts)
    for j=1:length(names)
        fprintf(of, '%s||6\t%s||6\t%i|\n',unique_clusts{i}, names{j}, sig_Assign(i,j)*1);
    end
    if i==1
       fprintf(of, '|- VertSpacer1|Spacer|3\n');
    end
end

fclose(of);



