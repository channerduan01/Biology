% Input IDX_MATRIX: 
%   rows number is repeat times,
%   columns number is items number of basic data
%
function ClusterAnalysis(target_cluster_idx, Q, protein_clusters, MRNA, PROTEIN_ORIGINAL)
gene_idxs = protein_clusters{target_cluster_idx};
figure();
hold on;
subplot(121), imagesc(MRNA(:,gene_idxs)');
subplot(122), imagesc(PROTEIN_ORIGINAL(:,gene_idxs)');
hold off;

related_cluster_idxs = zeros(1,length(gene_idxs));
fprintf('releated mrna clusters: ');
for i = 1:length(gene_idxs)
    [~,idx] = max(Q(:,gene_idxs(i)));
    related_cluster_idxs(i) = idx;
    fprintf('%d', idx);
    if i < length(gene_idxs)
        fprintf(', ');
    end
end
fprintf('\n');