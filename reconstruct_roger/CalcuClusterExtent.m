function [mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, Q, R_J, SELETION_THRESHOLD)
mrna_clusters = cell(K,1);
for i = 1:N
    [values, indices] = sort(Q(:,i));
    tmp_indices = indices(values > SELETION_THRESHOLD);
    for ii = 1:length(tmp_indices)
        mrna_clusters{tmp_indices(ii)} = [mrna_clusters{tmp_indices(ii)} i];
    end    
end
protein_clusters = cell(J,1);
for i = 1:N
    [values, indices] = sort(R_J(:,i));
    tmp_indices = indices(values > SELETION_THRESHOLD);
    for ii = 1:length(tmp_indices)
        protein_clusters{tmp_indices(ii)} = [protein_clusters{tmp_indices(ii)} i];
    end
end