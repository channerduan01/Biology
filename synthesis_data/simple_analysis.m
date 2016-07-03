%%
% Methods comparison

close all;

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

repeat_run = 2;

IDX_MATRIX_MRNA = zeros(2, N);
IDX_MATRIX_PROTEIN = zeros(2, N);

index_ = 1;
[~, idx] = max(H1_ORIGINAL);
IDX_MATRIX_MRNA(index_,:) = idx;
[~, idx] = max(H2_ORIGINAL);
IDX_MATRIX_PROTEIN(index_,:) = idx;


%% bench mark, k-means
tic
index_ = 2;
opts = statset('Display','final');
[idx1,~] = kmeans(MRNA',K,'Replicates',repeat_run,'Options',opts,'Start','sample');
IDX_MATRIX_MRNA(index_,:) = idx1;
H1 = zeros(K,length(idx1));
for i = 1:length(idx1)
    H1(idx1(i),i) = 1;
end
[idx2,~] = kmeans(PROTEIN',J,'Replicates',repeat_run,'Options',opts,'Start','sample');
IDX_MATRIX_PROTEIN(index_,:) = idx2;
H2 = zeros(J,length(idx2));
for i = 1:length(idx2)
    H2(idx2(i),i) = 1;
end
THETA1 = CalcuTheta(H1, H2, K, J, N);
% recover THETA to original order!
THETA1 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA1, K, J);
[~, entropy_j_k1, entropy_k_j1] = EntropyCalculate(K, J, sum(H1,2)/N, THETA1);
[mrna_correct,mrna_C1] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C1] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
% fprintf('k-means mrna entropy_j_k: %f, protein entropy_k_j: %f\n', entropy_j_k1, entropy_k_j1);
fprintf('k-means mrna correct rate: %f, protein correct rate: %f\n', mrna_correct, protein_correct);
% fprintf('k-means norm-fro of dist: %.2f, norm-1 of dist: %.2f\n', ...
%     norm(THETA1-THETA_ORIGINAL,'fro'),norm(THETA1-THETA_ORIGINAL,1));
toc

% drawCheckDataDistribution(MRNA',idx1,'MRNA K-MEANS');

% %% Roger's statistic model
% tic
% index_ = 3;
% [RESULT_, best_result, err_num] = ...
%     MyCoupleClusteringRepeat(MRNA, PROTEIN, K, J, 10, 10, true, repeat_run);
% [~, idx1] = max(best_result.Q);
% IDX_MATRIX_MRNA(index_,:) = idx1;
% [~, idx2] = max(best_result.Q_J);
% IDX_MATRIX_PROTEIN(index_,:) = idx2;
% THETA2 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, best_result.THETA, K, J);
% [~, entropy_j_k2, entropy_k_j2] = EntropyCalculate(K, J, sum(best_result.Q,2)/N, THETA2);
% [mrna_correct,mrna_C3] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
% [protein_correct,protein_C3] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
% % fprintf('Roger''s mrna entropy_j_k: %f, protein entropy_k_j: %f\n', entropy_j_k2, entropy_k_j2);
% fprintf('Roger''s mrna correct rate: %f, protein correct rate: %f\n', mrna_correct, protein_correct);
% % fprintf('Roger''s norm-fro of dist: %.2f, norm-1 of dist: %.2f\n', ...
% %     norm(THETA2-THETA_ORIGINAL,'fro'),norm(THETA2-THETA_ORIGINAL,1));
% toc
% 
% drawCheckDataDistribution(MRNA',idx1,'MRNA ROGER');


%% SNMF innovative model
tic
index_ = 4;
[W1,H1,W2,H2,THETA3,HIS,last_cost] = ...
    CoNMF(MRNA, PROTEIN, K, J, 0.1, 0.01, 0.01, 100, false, 0.001);
last_cost
[~, idx1] = max(H1);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(H2);
IDX_MATRIX_PROTEIN(index_,:) = idx2;
for k = 1:K
    THETA3(k,:) = THETA3(k,:)./sum(THETA3(k,:));
end
% recover THETA to original order!
THETA3 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA3, K, J);
[mrna_correct,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
fprintf('SNMF mrna correct rate: %f, protein correct rate: %f\n', mrna_correct, protein_correct);
% fprintf('SNMF norm-fro of dist: %.2f, norm-1 of dist: %.2f\n', ...
%     norm(THETA3-THETA_ORIGINAL,'fro'),norm(THETA3-THETA_ORIGINAL,1));
toc

% drawCheckDataDistribution(MRNA',idx1,'MRNA SNMF');


%% compare the patterns
figure();
set(gca,'FontSize',16);
hold on;
subplot(221), imagesc(THETA_ORIGINAL);
axis('off');
title('Original Theta', 'FontSize', 20)
subplot(222), imagesc(THETA1);
axis('off');
title('Coupled K-means Theta', 'FontSize', 20)
% subplot(223), imagesc(THETA2);
% axis('off');
% title('Roger''s Thetsa', 'FontSize', 20)
subplot(224), imagesc(THETA3);
axis('off');
title('Coupled SNMF Theta', 'FontSize', 20)
hold off;


% %% Output all errors of Theta
% res_ = zeros(1, repeat_run);
% for i = 1:repeat_run
%     [~, idx1] = max(RESULT_{i}.Q);
%     IDX_MATRIX_MRNA(index_,:) = idx1;
%     [~, idx2] = max(RESULT_{i}.Q_J);
%     IDX_MATRIX_PROTEIN(index_,:) = idx2;
%     THETA2 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), ...
%         idx2, RESULT_{i}.THETA, K, J);
%     fprintf('Roger''s norm-fro of dist: %.2f, norm-1 of dist: %.2f\n', ...
%         norm(THETA2-THETA_ORIGINAL,'fro'),norm(THETA2-THETA_ORIGINAL,1));
%     res_(i) = norm(THETA2-THETA_ORIGINAL,'fro');
% end
% 
% % best_result = RESULT_{48};






