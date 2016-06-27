%% Roger EM coupled clustering model
%
% init data
close all
clc
clear
[MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] = GeneDataLoad();
K = 15;
J = 19;



% %% K-mean test
% cost = @(V,W,H) norm(V-W*H,'fro');
% % protein = cell(J,1);
% % [centres, idx] = mykmeans(PROTEIN, J, 10);
% % for i = 1:N
% %     protein{idx(:,i)==1} = [protein{idx(:,i)==1} i];
% % end
% repeat = 1;
% opts = statset('Display','final');
% [idx,centres] = kmeans(PROTEIN',J,'Replicates',repeat,'Options',opts,'Start','sample');
% W = centres';
% H = zeros(J,length(idx));
% for i = 1:length(idx)
%     H(idx(i),i) = 1;
% end
% disp(['k-means cost:',num2str(cost(PROTEIN,W,H))]);

%% Roger's Model
MAX_ITER = 100;
patience = 1;
REPEAPT = 10;
RESULT = cell(REPEAPT,1);
index = 0;
error_result_num = 0;
while true
    % -------- permuted Proteins, keep annotation
    %     ii = randperm(N);
    %     MRNA(:,1:N) = MRNA(:,ii);
    %     PROTEIN(:,1:N) = PROTEIN(:,ii);
    % --------
    index = index + 1;
    if index > REPEAPT, break;end
    try
        [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = ...
            MyCoupleClustering(MRNA, PROTEIN, K, J, MAX_ITER, patience, true);
        [THETA_reverse, entropy_j_k, entropy_k_j] = EntropyCalculate(K, J, PI_K, THETA);
        [R_J, Q_J] = CalcuSubclusterBelonging(MRNA, AVG_K, VARIANCE_K, PROTEIN, AVG_J, VARIANCE_J, PI_K, THETA_reverse, R, K, J, N, T);
        low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N);
    catch ME
        low_bound = NaN;
        fprintf('error!!!\n')
    end
    if isnan(low_bound) || low_bound == -Inf || low_bound == Inf
        error_result_num = error_result_num + 1;
        index = index - 1;
        continue;
    end
    RESULT{index} = struct('low_bound',low_bound,'entropy_j_k',entropy_j_k,'entropy_k_j',entropy_k_j, ...
        'THETA_reverse',THETA_reverse,'Q',Q,'R',R,'R_J',R_J,'Q_J',Q_J,'PI_K',PI_K,'AVG_K',AVG_K, ...
        'VARIANCE_K',VARIANCE_K,'THETA',THETA,'AVG_J',AVG_J,'VARIANCE_J',VARIANCE_J);
    fprintf('Low bound is %f\n', low_bound);
    %     fprintf('cost is %f\n', cost(AVG_J,AVG_K,THETA_reverse'));
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
        [~, idx] = max(RESULT{i}.Q_J);
        IDX_MATRIX_PROTEIN(i,:) = idx;
    end
    [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
    [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
    fprintf('mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);
end
%% Single Analysis
SELETION_THRESHOLD = 0.3;
[mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, Q, Q_J, SELETION_THRESHOLD);
OutputClustersNM(K, J, mrna_clusters, protein_clusters, names);
ANALYSIS_PROTEIN_CLUSTER_IDX = 5;
ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, Q, protein_clusters, MRNA, PROTEIN_ORIGINAL);


%% test2
co_K_J = zeros(K, J, N);
for i = 1:N
    for k = 1:K
        for j = 1:J
            co_K_J(k,j,i) = Q(k,i)*Q_J(j,i);
        end
    end
end
K_J = sum(co_K_J,3)/N;
for j = 1:J
    K_J(:,j) = K_J(:,j)./PI_K(:);
end









