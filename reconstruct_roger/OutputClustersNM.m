function OutputClustersNM(K, J, mrna_clusters, protein_clusters, is_name, names)
fprintf('\nmRNA Clusters:\n')
for k = 1:K
    fprintf('%2d: ', k);
    tmp_indices = mrna_clusters{k};
    printGenesInCluster(tmp_indices, is_name, names);
    fprintf('\n');
end
fprintf('\nProtein Clusters:\n')
for j = 1:J
    fprintf('%2d: ', j);
    tmp_indices = protein_clusters{j};
    printGenesInCluster(tmp_indices, is_name, names);
    fprintf('\n');
end
end

function printGenesInCluster(tmp_indices, is_name, names)
for i = 1:length(tmp_indices)
    if is_name
        name_ = names.GeneName{tmp_indices(i)};
    else
        name_ = names.NM{tmp_indices(i)};
    end
    if ~strcmp(name_,'')
        fprintf('%s,', name_);
    end
end
end