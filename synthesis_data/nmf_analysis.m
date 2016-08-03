%%
% Methods comparison

close all;
addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

method = 'BP';
w_coef = 2;
h_coef = 2;
t_coef = 1;
max_iter = 100;
min_iter = 100;

[~, idx_mrna] = max(H1_ORIGINAL);
[~, idx_protein] = max(H2_ORIGINAL);

REPEAT_NUM = 20;
RESULT_TABLE = zeros(REPEAT_NUM, 4*3);
iter_used = zeros(REPEAT_NUM, 3);

IDX_MATRIX_MRNA = zeros(REPEAT_NUM, N, 4);
IDX_MATRIX_PROTEIN = zeros(REPEAT_NUM, N, 4);


%% bench mark, k-means
%     tic
for repeat_time = 1:REPEAT_NUM
    fprintf('\n k-means run on %d >>>>>>\n', repeat_time);
    repeat_run = 1;
    opts = statset('Display','final');
    [idx1,W1] = kmeans(MRNA',K,'Replicates',repeat_run,'Options',opts,'Start','sample');
    W1 = W1';
    H1 = zeros(K,length(idx1));
    for i = 1:length(idx1)
        H1(idx1(i),i) = 1;
    end
    [idx2,W2] = kmeans(PROTEIN',J,'Replicates',repeat_run,'Options',opts,'Start','sample');
    W2 = W2';
    H2 = zeros(J,length(idx2));
    for i = 1:length(idx2)
        H2(idx2(i),i) = 1;
    end
    IDX_MATRIX_MRNA(repeat_time,:,1) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,1) = idx2;
    THETA1 = CalcuTheta(H1, H2, K, J, N);
    % recover THETA to original order!
    THETA1 = ArrangeTheta(idx_mrna, idx1, idx_protein, idx2, THETA1, K, J);
    [~, entropy_j_k1, entropy_k_j1] = EntropyCalculate(K, J, sum(H1,2)/N, THETA1);
    mrna_correct = purity(idx1, idx_mrna);
    protein_correct = purity(idx2, idx_protein);
    theta_error = norm(THETA_ORIGINAL-THETA1,'fro');
    fprintf('k-means theta-error: %f, mrna correct rate: %f, protein correct rate: %f\n', ...
        theta_error, mrna_correct, protein_correct);
    RESULT_TABLE(repeat_time, 1:3) = [theta_error, mrna_correct, protein_correct];
end
%     toc

%% SNMF innovative model, two stage
%     tic
for repeat_time = 1:REPEAT_NUM
    fprintf('\n double NMF run on %d >>>>>>\n', repeat_time);
    index_ = 3;
    [W1,H1,W2,H2,~,HIS,last_iter] = ...
        CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
        , 'W_COEF', w_coef, 'H_COEF', h_coef, 'T_COEF', t_coef ...
        , 'VERBOSE', 0, 'METHOD', method ...
        , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
        );
    iter_used(repeat_time, 1) = iter_used(repeat_time, 1) + last_iter;
    [~, idx1] = max(H1);
    [~, idx2] = max(H2);
    IDX_MATRIX_MRNA(repeat_time,:,2) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,2) = idx2;
    THETA2 = CalcuTheta(H1, H2, K, J, N);
    % recover THETA to original order!
    THETA2 = ArrangeTheta(idx_mrna, idx1, idx_protein, idx2, THETA2, K, J);
    mrna_correct = purity(idx1, idx_mrna);
    protein_correct = purity(idx2, idx_protein);
    theta_error = norm(THETA_ORIGINAL-THETA2,'fro');
    fprintf('Naive NMFs theta-error: %f, mrna correct rate: %f, protein correct rate: %f\n', ...
        theta_error, mrna_correct, protein_correct);
    RESULT_TABLE(repeat_time, 4:6) = [theta_error, mrna_correct, protein_correct];
end
%     toc

%% Rogers's model
%     tic
[RESULT_, best_result, err_num] = ...
    MyCoupleClusteringRepeat(MRNA, PROTEIN, K, J, 100, 0.1, true, REPEAT_NUM);
for repeat_time = 1:REPEAT_NUM
    [~, idx1] = max(RESULT_{repeat_time}.Q);
    [~, idx2] = max(RESULT_{repeat_time}.Q_J);
    IDX_MATRIX_MRNA(repeat_time,:,3) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,3) = idx2;
    THETA3 = ArrangeTheta(idx_mrna, idx1, idx_protein, idx2, RESULT_{repeat_time}.THETA, K, J);
    mrna_correct = purity(idx1, idx_mrna);
    protein_correct = purity(idx2, idx_protein);
    theta_error = norm(THETA_ORIGINAL-THETA3,'fro');
    fprintf('Rogers''s theta-error: %f, mrna correct rate: %f, protein correct rate: %f, low_bound: %f\n', ...
        theta_error, mrna_correct, protein_correct, RESULT_{repeat_time}.low_bound);
    RESULT_TABLE(repeat_time, 7:9) = [theta_error, mrna_correct, protein_correct];
end
[~, idx1] = max(best_result.Q);
[~, idx2] = max(best_result.Q_J);
THETA3 = ArrangeTheta(idx_mrna, idx1, idx_protein, idx2, best_result.THETA, K, J);
%     toc

%% innovative model, My own flow NMF
%     tic
for repeat_time = 1:REPEAT_NUM
    fprintf('\n Coupled NMF run on %d >>>>>>\n', repeat_time);
    MRNA_ = MRNA;
    PROTEIN_ = PROTEIN;
    idx_mrna_ = idx_mrna;
    idx_protein_ = idx_protein;
    % -------- permuted Proteins, keep annotation
    %         ii = randperm(N);
    %         % only activate random on MRNA or PROTEIN
    %         % activate both means no permution!
    % %         MRNA_(:,1:N) = MRNA(:,ii);
    % %         idx_mrna_(:,1:N) = idx_mrna(:,ii);
    %         PROTEIN_(:,1:N) = PROTEIN(:,ii);
    %         idx_protein_(:,1:N) = idx_protein(:,ii);
    % --------
    [W1,H1,W2,H2,THETA4,HIS,last_iter] = ...
        CoNMF_v4_flow(MRNA_, PROTEIN_, K, J ...
        , 'W_COEF', w_coef, 'H_COEF', h_coef, 'T_COEF', t_coef ...
        , 'VERBOSE', 0, 'METHOD', method ...
        , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
        );
    [~, idx1] = max(H1);
    [~, idx2] = max(H2);
    IDX_MATRIX_MRNA(repeat_time,:,4) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,4) = idx2;
    iter_used(repeat_time, 3) = iter_used(repeat_time, 3) + last_iter;
    % recover THETA to original order!
    THETA4 = ArrangeTheta(idx_mrna_, idx1, idx_protein_, idx2, THETA4, K, J);
    mrna_correct = purity(idx1, idx_mrna_);
    protein_correct = purity(idx2, idx_protein_);
    theta_error = norm(THETA_ORIGINAL-THETA4,'fro');
    fprintf('FLow NMFs theta-error: %f, mrna correct rate: %f, protein correct rate: %f\n', ...
        theta_error, mrna_correct, protein_correct);
    RESULT_TABLE(repeat_time, 10:12) = [theta_error, mrna_correct, protein_correct];
end
%     toc

%% compare the patterns
figure();
set(gca,'FontSize',16);
hold on;
subplot(321), imagesc(THETA_ORIGINAL);
axis('off');
title('Original Theta', 'FontSize', 20)
subplot(322), imagesc(THETA1);
axis('off');
title('Double K-means Theta', 'FontSize', 20)
subplot(323), imagesc(THETA2);
axis('off');
title('Double NMFs Theta', 'FontSize', 20)
subplot(324), imagesc(THETA3);
axis('off');
title('Roger''s Theta', 'FontSize', 20)
subplot(325), imagesc(THETA4);
axis('off');
title('Coupled NMF Theta', 'FontSize', 20)
hold off;

%% output the performance
% choose different measurement
mean_err = mean(RESULT_TABLE);
std_err = std(RESULT_TABLE);
% mean_err = median(RESULT_TABLE);
% std_err = max(RESULT_TABLE) - mean_err;

fprintf('\n%d run\n', REPEAT_NUM);
mean(iter_used)
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,1));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,1));
fprintf('%-20s %5.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f-%.2f\n' ...
    , 'Double k-means:', mean_err(2)*100, std_err(2)*100 ...
    , mean_err(3)*100, std_err(3)*100 ...
    , mean_err(1), std_err(1), mrna_consistency, protein_consistency);
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,2));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,2));
fprintf('%-20s %5.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f-%.2f\n' ...
    , 'Double NMFs:', mean_err(5)*100, std_err(5)*100 ...
    , mean_err(6)*100, std_err(6)*100 ...
    , mean_err(4), std_err(4), mrna_consistency, protein_consistency);
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,3));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,3));
fprintf('%-20s %5.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f-%.2f\n' ...
    , 'Rogers''s:', mean_err(8)*100, std_err(8)*100 ...
    , mean_err(9)*100, std_err(9)*100 ...
    , mean_err(7), std_err(7), mrna_consistency, protein_consistency);
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,4));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,4));
fprintf('%-20s %5.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f-%.2f\n' ...
    , 'Coupled NMF: ', mean_err(11)*100, std_err(11)*100 ...
    , mean_err(12)*100, std_err(12)*100 ...
    , mean_err(10), std_err(10), mrna_consistency, protein_consistency);







