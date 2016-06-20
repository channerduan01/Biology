%% Roger EM coupled clustering model
%
% init data
close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

[MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] = GeneDataLoad();
K = 15;
J = 19;
MAX_ITER = 50;
REPEAT = 10;

%% taste
[W1,H1,W2,H2,THETA,HIS] = ...
    CoNMF(MRNA, PROTEIN, K, J, 0.1, 0.01, 0.01, MAX_ITER, true, 1);
plot(HIS(:,2:5));
set(gca,'FontSize',16);
xlabel('iteration', 'FontSize', 16);
ylabel('cost', 'FontSize', 16);
title('Gradient', 'FontSize', 20)
legend('W1 and H1', 'W2 and H2', 'W1 and THETA', 'W2 and inverse THETA');

%% Consistency
IDX_MATRIX_MRNA = zeros(REPEAT, N);
IDX_MATRIX_PROTEIN = zeros(REPEAT, N);
for i = 1:REPEAT
    [W1,H1,W2,H2,THETA,HIS] = ...
        CoNMF(MRNA, PROTEIN, K, J, 0.1, 0.01, 0.01, MAX_ITER, true, 1);
    [~, idx] = max(H1);
    IDX_MATRIX_MRNA(i,:) = idx;
    [~, idx] = max(H2);
    IDX_MATRIX_PROTEIN(i,:) = idx;
end
[mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
[protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
fprintf('mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);

%% Single Analysis
for i = 1:N
    H1(:,i) = H1(:,i)/sum(H1(:,i));
    H2(:,i) = H2(:,i)/sum(H2(:,i));
end

SELETION_THRESHOLD = 0.3;
[mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, H1, H2, SELETION_THRESHOLD);

OutputClustersNM(K, J, mrna_clusters, protein_clusters, names);
ANALYSIS_PROTEIN_CLUSTER_IDX = 1;
ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, H1, protein_clusters, MRNA, PROTEIN);






