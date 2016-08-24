%% NMF coupled clustering model
%
% init data
close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

load hmec;
MRNA = data{1};
P1 = data{2};
MRNA = MRNA';
PROTEIN = P1(:,1:7)';
PROTEIN_ORIGINAL = P1';
% different normalization methods, matters a lot!
MRNA = normalize(MRNA);
PROTEIN = normalize(PROTEIN);
% MRNA = (MRNA+1)/2;
% PROTEIN = (PROTEIN+1)/2;
% MRNA(MRNA<0) = 0;
% PROTEIN(PROTEIN<0) = 0;
% MRNA = abs(MRNA);
% PROTEIN = abs(PROTEIN);
% MRNA = pow2(MRNA);
% PROTEIN = pow2(PROTEIN);
% MRNA = normalize_v2(MRNA);
% PROTEIN = normalize_v2(PROTEIN);

[T,N] = size(MRNA);

K = 15;
J = 19;

%% Consistency
REPEAT = 10;
RESULT = cell(REPEAT,1);
IDX_MATRIX_MRNA = zeros(REPEAT, N);
IDX_MATRIX_PROTEIN = zeros(REPEAT, N);
err_num = 0;
i = 0;
HIS = zeros(REPEAT, 4);

% configuration
min_iter = 20;
max_iter = 100;

while true
%     %     -------- permuted Proteins, keep annotation
%         ii = randperm(N);
% %         MRNA(:,1:N) = MRNA(:,ii);
%         PROTEIN(:,1:N) = PROTEIN(:,ii);
% %     --------
    
    i = i + 1;
    if i > REPEAT, break; end
    fprintf('process-iter >>>>>> %d\n', i);
    %     try
    %         [W1,H1,W2,H2,THETA,HIS] = ...
    %             CoNMF(MRNA, PROTEIN, 12, 12, 0.1, 0.01, 0.01, 100, true, 0.01);
    %     [W1_res,H1_res,W2_res,H2_res,~,HIS,last_iter] = CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
    %         , 'MAX_ITER', max_iter, 'MIN_ITER', max_iter, 'VERBOSE', 1, 'METHOD', 'BP' ...
    %         , 'W_COEF', wCoef, 'H_COEF', hCoef, 'T_COEF', 1, 'PATIENCE', 0.001 ...
    %         );
    %     Theta = CalcuTheta(H1_res, H2_res, K, J, N);
    [W1_res,H1_res,W2_res,H2_res,Theta,HIS,last_iter] = CoNMF_v4_flow(MRNA, PROTEIN, K, J ...
        , 'MAX_ITER', max_iter, 'MIN_ITER', min_iter, 'VERBOSE', 1, 'METHOD', 'BP' ...
        , 'W_COEF', 1, 'H_COEF', 1, 'T_COEF', 1, 'PATIENCE', 10 ...
        );
    PI_K = sum(normalizeColumn(H1_res),2)/N;
    [THETA_reverse, entropy_j_k, entropy_k_j] = EntropyCalculate(K, J, PI_K, Theta);
    entropy_j_k
    if sum(sum(isnan(Theta))) || isnan(entropy_j_k) > 0
        err_num = err_num+1;
        i = i - 1;
        continue;
    end
    %     catch
    %         err_num = err_num+1;
    %         i = i - 1;
    %         continue;
    %     end
    
    [~, idx] = max(H1_res);
    IDX_MATRIX_MRNA(i,:) = idx;
    [~, idx] = max(H2_res);
    IDX_MATRIX_PROTEIN(i,:) = idx;

    RESULT{i} = struct('W1',W1_res,'H1',H1_res,'W2',W2_res,'H2',H2_res, ...
        'entropy_j_k',entropy_j_k,'entropy_k_j',entropy_k_j, ...
        'Theta',Theta,'HIS',HIS,'last_iter',last_iter);
end

if REPEAT > 1
    [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
    [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
    fprintf('error: %d, mrna consistency: %f, protein consistency: %f\n', ...
        err_num, mrna_consistency, protein_consistency);
end

% %% trans
% A = W1_res;
% for i = 1:size(A,1)
%     A(i,:) = A(i,:)/sum(A(i,:));
% end
% H1_res = A';
% A = W2_res;
% for i = 1:size(A,1)
%     A(i,:) = A(i,:)/sum(A(i,:));
% end
% H2_res = A';
% THETA = CalcuTheta(H1_res, H2_res, K, J, N);
%% Draw the process of clustering
% plot(HIS(:,:));
% set(gca,'FontSize',16);
% xlabel('iteration', 'FontSize', 16);
% ylabel('cost', 'FontSize', 16);
% title('Gradient', 'FontSize', 20)
% legend('W1 and H1', 'W2 and H2', 'W1 and THETA', 'W2 and inverse THETA');

%% Single Analysis
for i = 1:N
    H1_res(:,i) = H1_res(:,i)/sum(H1_res(:,i));
    H2_res(:,i) = H2_res(:,i)/sum(H2_res(:,i));
end

SELETION_THRESHOLD = 0.5;
% [mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, H1_res, H2_res, SELETION_THRESHOLD);
[~, idx] = max(H1_res);
mrna_clusters = cell(K,1);
for k = 1:K
    mrna_clusters{k} = find(idx==k);
end
[~, idx] = max(H2_res);
protein_clusters = cell(J,1);
for j = 1:J
    protein_clusters{j} = find(idx==j);
end

OutputClustersNM(K, J, mrna_clusters, protein_clusters, true, names);
ANALYSIS_PROTEIN_CLUSTER_IDX = 5;
ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, H1_res, protein_clusters, MRNA, PROTEIN, names);

% %% Show the distribution of data source
% figure();
% hold on;
% subplot(121);
% hist(MRNA(:));
% set(gca,'FontSize',16);
% xlabel('express level', 'FontSize', 16);
% ylabel('quantity', 'FontSize', 16);
% title('Original MRNA data', 'FontSize', 20);
%
% subplot(122);
% hist(PROTEIN(:));
% set(gca,'FontSize',16);
% xlabel('express level', 'FontSize', 16);
% ylabel('quantity', 'FontSize', 16);
% title('Original protein data', 'FontSize', 20);
% hold off;


%% converge checking
ALL_HIS = zeros(max_iter+1,REPEAT);
for i = 1:REPEAT
    ALL_HIS(:,i) = sum(RESULT{i}.HIS,2);
end
figure();
plot(sum(ALL_HIS,2));
xlabel('iterations', 'FontSize', 16);
ylabel('cost', 'FontSize', 16);
title('Converge process', 'FontSize', 20)




