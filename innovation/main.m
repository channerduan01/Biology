%% NMF coupled clustering model
%
% init data
close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

[MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] = GeneDataLoad();
K = 15;
J = 19;

%% taste
[W1,H1,W2,H2,THETA,HIS] = ...
    CoNMF(MRNA, PROTEIN, K, J, 0.1, 0.01, 0.01, 100, true, 0.1);
plot(HIS(:,2:5));
set(gca,'FontSize',16);
xlabel('iteration', 'FontSize', 16);
ylabel('cost', 'FontSize', 16);
title('Gradient', 'FontSize', 20)
legend('W1 and H1', 'W2 and H2', 'W1 and THETA', 'W2 and inverse THETA');

for k = 1:K
    THETA(k,:) = THETA(k,:)/sum(THETA(k,:));
end

%% Consistency
MAX_ITER = 100;
REPEAT = 10;
patience = 1;

IDX_MATRIX_MRNA = zeros(REPEAT, N);
IDX_MATRIX_PROTEIN = zeros(REPEAT, N);
err_num = 0;
i = 0;
while true
    i = i + 1;
    if i > REPEAT, break; end
    try
    [W1,H1,W2,H2,THETA,HIS] = ...
        CoNMF(MRNA, PROTEIN, K, J, 1, 1, 1, MAX_ITER, true, patience);
    catch
        err_num = err_num+1;
        i = i - 1;
        continue;
    end
    [~, idx] = max(H1);
    IDX_MATRIX_MRNA(i,:) = idx;
    [~, idx] = max(H2);
    IDX_MATRIX_PROTEIN(i,:) = idx;
end
[mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
[protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
fprintf('error: %d, mrna consistency: %f, protein consistency: %f\n', ...
    err_num, mrna_consistency, protein_consistency);

% %% Single Analysis
% for i = 1:N
%     H1(:,i) = H1(:,i)/sum(H1(:,i));
%     H2(:,i) = H2(:,i)/sum(H2(:,i));
% end
% 
% SELETION_THRESHOLD = 0.5;
% [mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, H1, H2, SELETION_THRESHOLD);
% 
% OutputClustersNM(K, J, mrna_clusters, protein_clusters, names);
% ANALYSIS_PROTEIN_CLUSTER_IDX = 19;
% ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, H1, protein_clusters, MRNA, PROTEIN_ORIGINAL);






