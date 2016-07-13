%%
% Methods comparison

close all;

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

repeat_run = 1;

method = 'MU';


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

%% SNMF innovative model, two stage
tic
index_ = 3;
[W1,H1,W2,H2,~,HIS,last_iter] = ...
    CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
    , 'W_COEF', 2, 'H_COEF', 2, 'T_COEF', 1 ...
    , 'VERBOSE', 0, 'METHOD', method ...
    );
THETA2 = CalcuTheta(H1, H2, K, J, N);
[~, idx1] = max(H1);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(H2);
IDX_MATRIX_PROTEIN(index_,:) = idx2;
for k = 1:K
    THETA2(k,:) = THETA2(k,:)./sum(THETA2(k,:));
end
% recover THETA to original order!
THETA2 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA2, K, J);
[mrna_correct,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
fprintf('2 NMFs mrna correct rate: %f, protein correct rate: %f\n', mrna_correct, protein_correct);
% fprintf('SNMF norm-fro of dist: %.2f, norm-1 of dist: %.2f\n', ...
%     norm(THETA3-THETA_ORIGINAL,'fro'),norm(THETA3-THETA_ORIGINAL,1));
toc

% drawCheckDataDistribution(MRNA',idx1,'MRNA SNMF');

%% SNMF innovative model, coupled model
tic
index_ = 4;
[W1,H1,W2,H2,THETA3,HIS,last_iter] = ...
    CoNMF_v3_co2(MRNA, PROTEIN, K, J, 'MAX_ITER', 100, 'PATIENCE', 1 ...
    , 'W_COEF', 2, 'H_COEF', 2, 'T_COEF', 1 ...
    , 'VERBOSE', 0, 'METHOD', method ...
    );
[~, idx1] = max(H1);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(H2);
IDX_MATRIX_PROTEIN(index_,:) = idx2;
% recover THETA to original order!
THETA3 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA3, K, J);
[mrna_correct,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
fprintf('Coupled NMF mrna correct rate: %f, protein correct rate: %f\n', mrna_correct, protein_correct);
% fprintf('SNMF norm-fro of dist: %.2f, norm-1 of dist: %.2f\n', ...
%     norm(THETA3-THETA_ORIGINAL,'fro'),norm(THETA3-THETA_ORIGINAL,1));
toc

%% SNMF innovative model, flow model
tic
index_ = 5;
[W1,H1,W2,H2,THETA4,HIS,last_iter] = ...
    CoNMF_v4_flow(MRNA, PROTEIN, K, J ...
    , 'W_COEF', 2, 'H_COEF', 2, 'T_COEF', 1 ...
    , 'VERBOSE', 0, 'METHOD', method ...
    );
[~, idx1] = max(H1);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(H2);
IDX_MATRIX_PROTEIN(index_,:) = idx2;
% recover THETA to original order!
THETA4 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA4, K, J);
[mrna_correct,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
fprintf('Flow NMF mrna correct rate: %f, protein correct rate: %f\n', mrna_correct, protein_correct);
% fprintf('SNMF norm-fro of dist: %.2f, norm-1 of dist: %.2f\n', ...
%     norm(THETA3-THETA_ORIGINAL,'fro'),norm(THETA3-THETA_ORIGINAL,1));
toc


%% compare the patterns
figure();
set(gca,'FontSize',16);
hold on;
subplot(321), imagesc(THETA_ORIGINAL);
axis('off');
title('Original Theta', 'FontSize', 20)
subplot(322), imagesc(THETA1);
axis('off');
title('Coupled K-means Theta', 'FontSize', 20)
subplot(323), imagesc(THETA2);
axis('off');
% title('Roger''s Theta', 'FontSize', 20)
title('Naive two NMFs Theta', 'FontSize', 20)
subplot(324), imagesc(THETA3);
axis('off');
title('Coupled NMF Theta', 'FontSize', 20)
subplot(325), imagesc(THETA4);
axis('off');
title('Flow NMF Theta', 'FontSize', 20)
hold off;





