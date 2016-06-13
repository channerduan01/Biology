% Roger EM coupled clustering model
%
% init data
close all
clear
clc

cost = @(V,W,H) norm(V-W*H,'fro');

load hmec;
MRNA = data{1};
P1 = data{2}; 
MRNA = normalize(MRNA)';
PROTEIN = normalize(P1(:,2:7))';
PROTEIN_ORIGINAL = normalize(P1)';
[T,N] = size(MRNA);
K = 15;
J = 19;

%% K-mean test
% protein = cell(J,1);
% [centres, idx] = mykmeans(PROTEIN, 19, 10);
% for i = 1:N
%     protein{idx(:,i)==1} = [protein{idx(:,i)==1} i];
% end
% % repeat = 1;
% % opts = statset('Display','final');
% % [idx,centres] = kmeans(V',K,'Replicates',repeat,'Options',opts,'Start','sample');
% % W = centres';
% % H = zeros(K,length(idx));
% % for i = 1:length(idx)
% %     H(idx(i),i) = 1;
% % end
% % disp(['k-means cost:',num2str(cost(V,W,H))]);

%% Roger's Model
MAX_ITER = 100;
patience = 2;
REPEAPT = 1;
RESULT = cell(REPEAPT,1);
for i = 1:REPEAPT
    [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = ...
        MyCoupleClustering(MRNA, PROTEIN, K, J, MAX_ITER, patience, true);
    [THETA_reverse, entropy_j_k, entropy_k_j] = BasicCoupleClusteringAnalysis(K, J, PI_K, THETA);
    low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N);
    R_J = mean(R,3);
    RESULT{i} = struct('low_bound',low_bound,'entropy_j_k',entropy_j_k,'entropy_k_j',entropy_k_j, ...
        'THETA_reverse',THETA_reverse,'Q',Q,'R',R,'R_J',R_J,'PI_K',PI_K,'AVG_K',AVG_K, ...
        'VARIANCE_K',VARIANCE_K,'THETA',THETA,'AVG_J',AVG_J,'VARIANCE_J',VARIANCE_J);
    fprintf('Low bound is %f\n', low_bound);
    fprintf('Entropy of p(j|k) = %f\n', entropy_j_k);
    fprintf('Entropy of p(k|j) = %f\n', entropy_k_j);
end
%% Consistency
if REPEAPT > 1
    IDX_MATRIX_MRNA = zeros(REPEAPT, N);
    IDX_MATRIX_PROTEIN = zeros(REPEAPT, N);
    for i = 1:REPEAPT
        [~, idx] = max(RESULT{i}.Q);
        IDX_MATRIX_MRNA(i,:) = idx;
        [~, idx] = max(RESULT{i}.R_J);
        IDX_MATRIX_PROTEIN(i,:) = idx;  
    end
    [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
    [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
    fprintf('mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);
end
%% Single Analysis
SELETION_THRESHOLD = 0.4;
[mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, Q, R_J, SELETION_THRESHOLD);
OutputClustersNM(K, J, mrna_clusters, protein_clusters, names);
ANALYSIS_PROTEIN_CLUSTER_IDX = 12;
ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, Q, protein_clusters, MRNA, PROTEIN);

    
    
    
    


