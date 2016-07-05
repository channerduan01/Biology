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
PROTEIN_ORIGINAL = normalize(P1)';
[T,N] = size(MRNA);

K = 15;
J = 19;

% %% taste
% [W1,H1,W2,H2,THETA,HIS] = ...
%     CoNMF(MRNA, PROTEIN, K, J, 0.1, 0.01, 0.01, 100, true, 0.01);
% % plot(HIS(:,2:5));
% % set(gca,'FontSize',16);
% % xlabel('iteration', 'FontSize', 16);
% % ylabel('cost', 'FontSize', 16);
% % title('Gradient', 'FontSize', 20)
% % legend('W1 and H1', 'W2 and H2', 'W1 and THETA', 'W2 and inverse THETA');
%
% for k = 1:K
%     THETA(k,:) = THETA(k,:)/sum(THETA(k,:));
% end

%% Consistency
REPEAT = 2;

IDX_MATRIX_MRNA = zeros(REPEAT, N);
IDX_MATRIX_PROTEIN = zeros(REPEAT, N);
err_num = 0;
i = 0;
while true
    i = i + 1;
    if i > REPEAT, break; end
    fprintf('process-iter >>>>>> %d\n', i);
    %     try
    %         [W1,H1,W2,H2,THETA,HIS] = ...
    %             CoNMF(MRNA, PROTEIN, 12, 12, 0.1, 0.01, 0.01, 100, true, 0.01);
    
    [W1,H1,W2,H2,~,HIS,~] = ...
        CoNMF_v2_separate(MRNA, PROTEIN, K, J, 'MAX_ITER', 50 ...
        , 'W_COEF', 2, 'H_COEF', 2, 'T_COEF', 2 ...
        ...%    , 'W1_INIT', W1, 'H1_INIT', H1 ...
        ...%     , 'W2_INIT', W2, 'H2_INIT', H2 ...
        );
    
    %     [W1,H1,W2,H2,~,~,~] = ...
    %         CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
    %         , 'W_COEF', 2, 'H_COEF', 2, 'T_COEF', 1 ...
    %         , 'MAX_ITER', 50, 'verbose', 1 ...
    %         ...%    , 'W1_INIT', W1, 'H1_INIT', H1 ...
    %         ...%     , 'W2_INIT', W2, 'H2_INIT', H2 ...
    %         );
    
    THETA = CalcuTheta(H1, H2, K, J, N);
    
    %     catch
    %         err_num = err_num+1;
    %         i = i - 1;
    %         continue;
    %     end
    [~, idx] = max(H1);
    IDX_MATRIX_MRNA(i,:) = idx;
    [~, idx] = max(H2);
    IDX_MATRIX_PROTEIN(i,:) = idx;
end
if REPEAT > 1
    [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
    [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
    fprintf('error: %d, mrna consistency: %f, protein consistency: %f\n', ...
        err_num, mrna_consistency, protein_consistency);
end

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






