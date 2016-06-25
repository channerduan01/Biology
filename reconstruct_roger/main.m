%% Roger EM coupled clustering model
%
% init data
close all
clear
clc

[MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] = GeneDataLoad();
K = 15;
J = 19;

cost = @(V,W,H) norm(V-W*H,'fro');

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
patience = 10;
REPEAPT = 1;
RESULT = cell(REPEAPT,1);
i = 0;
error_result_num = 0;
while true
    % -------- permuted Proteins, keep annotation
%     ii = randperm(N);
%     MRNA(:,1:N) = MRNA(:,ii);
%     PROTEIN(:,1:N) = PROTEIN(:,ii);
    % --------
    i = i + 1;
    if i > REPEAPT, break;end
    [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = ...
        MyCoupleClustering(MRNA, PROTEIN, K, J, MAX_ITER, patience, true);
    R_J = zeros(J, N);
    for j = 1:J
        for i = 1:N
            tmp = R(j,i,:);
            R_J(j,i) = sum(tmp(:).*PI_K);
        end
    end
    [THETA_reverse, entropy_j_k, entropy_k_j] = EntropyCalculate(K, J, PI_K, THETA);
    low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N);
    if isnan(low_bound) || low_bound == -Inf || low_bound == Inf
        error_result_num = error_result_num + 1;
        i = i - 1;
        continue;
    end
    RESULT{i} = struct('low_bound',low_bound,'entropy_j_k',entropy_j_k,'entropy_k_j',entropy_k_j, ...
        'THETA_reverse',THETA_reverse,'Q',Q,'R',R,'R_J',R_J,'PI_K',PI_K,'AVG_K',AVG_K, ...
        'VARIANCE_K',VARIANCE_K,'THETA',THETA,'AVG_J',AVG_J,'VARIANCE_J',VARIANCE_J);
    fprintf('Low bound is %f\n', low_bound);
    fprintf('cost is %f\n', cost(AVG_J,AVG_K,THETA_reverse'));
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
SELETION_THRESHOLD = 0.3;
[mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, Q, Q_J, SELETION_THRESHOLD);
OutputClustersNM(K, J, mrna_clusters, protein_clusters, names);
ANALYSIS_PROTEIN_CLUSTER_IDX = 5;
ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, Q, protein_clusters, MRNA, PROTEIN_ORIGINAL);
 

%% test1
PI_J = sum(R_J,2)/N;
Q_J = zeros(J, N);
for j = 1:J
   for i = 1:N
       for k = 1:K
        Q_J(j,i) = Q_J(j,i) + PI_J(j)*THETA_reverse(j,k)*mvnpdf(MRNA(:,i)',AVG_K(:,k)',VARIANCE_K(k)*eye(T))...
            *mvnpdf(PROTEIN(:,i)',AVG_J(:,j)',VARIANCE_J(j)*eye(T));
       end
   end
end
for i = 1:N
    Q_J(:,i) = Q_J(:,i)/sum(Q_J(:,i));
end

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









