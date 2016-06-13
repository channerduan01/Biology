function OutputClustersNM(K, J, mrna_clusters, protein_clusters, names)
fprintf('\nmRNA Clusters:\n')
for k = 1:K
    fprintf('%d: ', k);
    tmp_indices = mrna_clusters{k};
    for ii = 1:length(tmp_indices)
        fprintf('%s,', names.GeneName{tmp_indices(ii)});
    end
    fprintf('\n');
end
fprintf('\nProtein Clusters:\n')
for j = 1:J
    fprintf('%d: ', j);
    tmp_indices = protein_clusters{j};
    for ii = 1:length(tmp_indices)
        fprintf('%s,', names.GeneName{tmp_indices(ii)});
    end
    fprintf('\n');
end