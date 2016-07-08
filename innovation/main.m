%% NMF coupled clustering model
%
% init data
close all
% clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

load hmec;
MRNA = data{1};
P1 = data{2};
MRNA = MRNA';
PROTEIN = P1(:,2:7)';
PROTEIN_ORIGINAL = P1';
% different normalization methods, matters a lot!
% MRNA = normalize(MRNA);
% PROTEIN = normalize(PROTEIN);
% MRNA(MRNA<0) = 0;
% PROTEIN(PROTEIN<0) = 0;
MRNA = pow2(MRNA);
PROTEIN = pow2(PROTEIN);
% =========================== how to proper normalize ?

[T,N] = size(MRNA);

K = 15;
J = 19;


%% Consistency
REPEAT = 100;

IDX_MATRIX_MRNA = zeros(REPEAT, N);
IDX_MATRIX_PROTEIN = zeros(REPEAT, N);

err_num = 0;
i = 0;
HIS = zeros(REPEAT, 4);

% configuration
wCoef = 1;
hCoef = 1;
max_iter = 100;

while true
    i = i + 1;
    if i > REPEAT, break; end
    fprintf('process-iter >>>>>> %d\n', i);
    %     try
    %         [W1,H1,W2,H2,THETA,HIS] = ...
    %             CoNMF(MRNA, PROTEIN, 12, 12, 0.1, 0.01, 0.01, 100, true, 0.01);
    [W1_res,H1_res,W2_res,H2_res,~,HIS,last_iter] = CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
        , 'MAX_ITER', max_iter, 'MIN_ITER', max_iter, 'verbose', 1 ...
        , 'W_COEF', wCoef, 'H_COEF', hCoef, 'T_COEF', 1 ...
        );
    
    THETA = CalcuTheta(H1_res, H2_res, K, J, N);
    if sum(sum(isnan(THETA))) > 0
        err_num = err_num+1;
        i = i - 1;
        continue;
    end
    for k = 1:K
        THETA(k,:) = THETA(k,:)/sum(THETA(k,:));
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
end

if REPEAT > 1
    [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
    [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
    fprintf('error: %d, mrna consistency: %f, protein consistency: %f\n', ...
        err_num, mrna_consistency, protein_consistency);
end

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

SELETION_THRESHOLD = 0.3;
[mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, H1_res, H2_res, SELETION_THRESHOLD);

OutputClustersNM(K, J, mrna_clusters, protein_clusters, names);
ANALYSIS_PROTEIN_CLUSTER_IDX = 10;
ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, H1_res, protein_clusters, MRNA, PROTEIN);






