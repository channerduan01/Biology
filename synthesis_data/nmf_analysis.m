%%
% Methods comparison

close all;
addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));


max_iter = 100;
min_iter = 20;

[~, idx_mrna] = max(H1_ORIGINAL);
[~, idx_protein] = max(H2_ORIGINAL);

REPEAT_NUM = 20;
RESULT_TABLE = zeros(REPEAT_NUM, 4*5);
iter_used = zeros(REPEAT_NUM, 3);

IDX_MATRIX_MRNA = zeros(REPEAT_NUM, N, 4);
IDX_MATRIX_PROTEIN = zeros(REPEAT_NUM, N, 4);

%%
method = 'BP';
w_coef = 6;
h_coef = 6;

%% benchmark, k-means
%     tic
RESULT_KMEAN = cell(REPEAT_NUM,1);
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
    total_cost = norm(MRNA-W1*H1,'fro') + norm(PROTEIN-W2*H2,'fro');
    IDX_MATRIX_MRNA(repeat_time,:,1) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,1) = idx2;
    THETA1 = CalcuTheta(H1, H2, K, J, N);
    % recover THETA to original order!
    THETA1 = ArrangeTheta(idx_mrna, idx1, idx_protein, idx2, THETA1, K, J);
    [~, entropy_j_k1, entropy_k_j1] = EntropyCalculate(K, J, sum(H1,2)/N, THETA1);
    mrna_correct = purity(idx1, idx_mrna);
    protein_correct = purity(idx2, idx_protein);
    theta_error = norm(THETA_ORIGINAL-THETA1,'fro');
    fprintf('k-means theta-error: %f, mrna correct rate: %f, protein correct rate: %f, total_cost: %f\n', ...
        theta_error, mrna_correct, protein_correct, total_cost);
    RESULT_KMEAN{repeat_time} = struct('total_cost',total_cost,'mrna_correct',mrna_correct, ...
        'protein_correct', protein_correct, 'theta_error',theta_error, 'theta', THETA1);
    RESULT_TABLE(repeat_time, 1:5) = [theta_error, mrna_correct, protein_correct, mean(sparsity(H1)), mean(sparsity(H2))];
end
best_kmeans_idx = 1;
for repeat_time = 2:REPEAT_NUM
    if RESULT_KMEAN{repeat_time}.total_cost < RESULT_KMEAN{best_kmeans_idx}.total_cost
        best_kmeans_idx = repeat_time;
    end
end
THETA1 = RESULT_KMEAN{best_kmeans_idx}.theta;

%     toc

%% NMF innovative model, two stage
%     tic
RESULT_DNMF = cell(REPEAT_NUM,1);
repeat_time = 1;
while true
    fprintf('\n double NMF run on %d >>>>>>\n', repeat_time);
    index_ = 3;
%     [W1,H1,W2,H2,~,HIS,last_iter] = ...
%         CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
%         , 'W_COEF', w_coef, 'H_COEF', h_coef, 'T_COEF', t_coef ...
%         , 'VERBOSE', 0, 'METHOD', method ...
%         , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
%         );
    [W1,H1,W2,H2,~,HIS,last_iter] = ...
        CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
        , 'W_COEF', w_coef, 'H_COEF', h_coef ...
        , 'VERBOSE', 0, 'METHOD', method, 'PATIENCE', 0.1 ...
        , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
        );    
    total_cost = sum(HIS(last_iter+1,1:2));
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
    if isnan(total_cost) || sum(isnan(theta_error)) > 0
       continue; 
    end
    fprintf('Naive NMFs theta-error: %f, mrna correct rate: %f, protein correct rate: %f, total_cost: %f\n', ...
        theta_error, mrna_correct, protein_correct, total_cost);
    RESULT_DNMF{repeat_time} = struct('total_cost',total_cost,'mrna_correct',mrna_correct, ...
        'protein_correct', protein_correct, 'theta_error',theta_error, 'theta', THETA2);
    RESULT_TABLE(repeat_time, 6:10) = [theta_error, mrna_correct, protein_correct, mean(sparsity(H1)), mean(sparsity(H2))];
    repeat_time = repeat_time + 1;
    if repeat_time > REPEAT_NUM, break; end
end
best_dnmf_idx = 1;
for repeat_time = 2:REPEAT_NUM
    if RESULT_DNMF{repeat_time}.total_cost < RESULT_DNMF{best_dnmf_idx}.total_cost
        best_dnmf_idx = repeat_time;
    end
end
THETA2 = RESULT_DNMF{best_dnmf_idx}.theta;
%     toc

%% Rogers's model
% THETA3 = THETA2; % just for test
%     tic
[RESULT_, best_result, err_num] = ...
    MyCoupleClusteringRepeat(MRNA, PROTEIN, K, J, 100, 1, true, REPEAT_NUM);
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
    fprintf('%f\n', mrna_correct+protein_correct);
    RESULT_TABLE(repeat_time, 11:15) = [theta_error, mrna_correct, protein_correct, ...
        mean(sparsity(RESULT_{repeat_time}.Q)), mean(sparsity(RESULT_{repeat_time}.Q_J))];
end
[~, idx1] = max(best_result.Q);
[~, idx2] = max(best_result.Q_J);
THETA3 = ArrangeTheta(idx_mrna, idx1, idx_protein, idx2, best_result.THETA, K, J);
best_rogers_idx = 1;
for repeat_time = 1:REPEAT_NUM
    if RESULT_{repeat_time}.low_bound == best_result.low_bound
        best_rogers_idx = repeat_time;
        break;
    end
end

%     toc

%% innovative model, My own coupled NMF
%     tic
RESULT_CNMF = cell(REPEAT_NUM,1);
repeat_time = 1;
while true
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
    [W1,H1,W2,H2,THETA4,HIS,last_iter,cost_nmf] = ...
        CoNMF_v4_flow(MRNA_, PROTEIN_, K, J ...
        , 'W_COEF', w_coef, 'H_COEF', h_coef ...
        , 'T_COEF', 1 ...
        , 'VERBOSE', 0, 'METHOD', method, 'PATIENCE', 1 ...
        , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
        );
    [~, idx1] = max(H1);
    [~, idx2] = max(H2);
    IDX_MATRIX_MRNA(repeat_time,:,4) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,4) = idx2;
    total_cost = sum(HIS(last_iter+1,1:2));
    iter_used(repeat_time, 3) = iter_used(repeat_time, 3) + last_iter;
    % recover THETA to original order!
    THETA4 = ArrangeTheta(idx_mrna_, idx1, idx_protein_, idx2, THETA4, K, J);
    mrna_correct = purity(idx1, idx_mrna_);
    protein_correct = purity(idx2, idx_protein_);
    theta_error = norm(THETA_ORIGINAL-THETA4,'fro');
    if isnan(total_cost) || sum(isnan(theta_error)) > 0
       continue; 
    end
    fprintf('Coupled NMFs cost: %f, theta-error: %f, mrna correct rate: %f, protein correct rate: %f, total_cost: %f\n', ...
        cost_nmf, theta_error, mrna_correct, protein_correct, total_cost);
    RESULT_CNMF{repeat_time} = struct('total_cost',total_cost,'mrna_correct',mrna_correct, ...
        'protein_correct', protein_correct, 'theta_error',theta_error, 'theta', THETA4);
    RESULT_TABLE(repeat_time, 16:20) = [theta_error, mrna_correct, protein_correct, mean(sparsity(H1)), mean(sparsity(H2))];
    repeat_time = repeat_time + 1;
    if repeat_time > REPEAT_NUM, break; end
end

best_cnmf_idx = 1;
for repeat_time = 2:REPEAT_NUM
    if RESULT_CNMF{repeat_time}.total_cost < RESULT_CNMF{best_cnmf_idx}.total_cost
        best_cnmf_idx = repeat_time;
    end
end
THETA4 = RESULT_CNMF{best_cnmf_idx}.theta;
%     toc


% %% real coupled k-means
%     [W1,H1,W2,H2,THETA6,cost_] = ...
%         CoupledKmeans(MRNA_, PROTEIN_, K, J, 0.9, 30);

%% compare the patterns
figure();
set(gca,'FontSize',16);
hold on;
subplot(151), imagesc(THETA_ORIGINAL);
axis('off');
title('Real Correlations', 'FontSize', 21)
subplot(152), imagesc(THETA1);
axis('off');
title('Double K-means', 'FontSize', 21)
subplot(153), imagesc(THETA2);
axis('off');
title('Double NMF', 'FontSize', 21)
subplot(154), imagesc(THETA3);
axis('off');
title('Rogers'' model', 'FontSize', 21)
subplot(155), imagesc(THETA4);
axis('off');
title('Coupled NMF', 'FontSize', 21)
hold off;

%% output the performance
% choose different measurement
mean_err = mean(RESULT_TABLE);
std_err = std(RESULT_TABLE);

% mean_err = median(RESULT_TABLE);
% std_err = max(RESULT_TABLE) - mean_err;

fprintf('\n%d run\n', REPEAT_NUM);
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,1));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,1));
fprintf('%-20s , %.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f + %.2f, %.0f%%, %.0f%%, %.2f + %.2f\n' ...
    , 'Double K-means', mean_err(2)*100, std_err(2)*100 ...
    , mean_err(3)*100, std_err(3)*100 ...
    , mean_err(4), std_err(4) ...
    , mean_err(5), std_err(5) ...
    , mrna_consistency*100, protein_consistency*100 ...
    , mean_err(1), std_err(1));
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,2));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,2));
fprintf('%-20s , %.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f + %.2f, %.0f%%, %.0f%%, %.2f + %.2f\n' ...
    , 'Double NMFs', mean_err(7)*100, std_err(7)*100 ...
    , mean_err(8)*100, std_err(8)*100 ...
    , mean_err(9), std_err(9) ...
    , mean_err(10), std_err(10) ...
    , mrna_consistency*100, protein_consistency*100 ...
    , mean_err(6), std_err(6));
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,3));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,3));
fprintf('%-20s , %.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f + %.2f, %.0f%%, %.0f%%, %.2f + %.2f\n' ...
    , 'Rogers'' Model', mean_err(12)*100, std_err(12)*100 ...
    , mean_err(13)*100, std_err(13)*100 ...
    , mean_err(14), std_err(14) ...
    , mean_err(15), std_err(15) ...
    , mrna_consistency*100, protein_consistency*100 ...
    , mean_err(11), std_err(11));
[mrna_consistency,~] = CalcuConsistency(IDX_MATRIX_MRNA(:,:,4));
[protein_consistency,~] = CalcuConsistency(IDX_MATRIX_PROTEIN(:,:,4));
fprintf('%-20s , %.1f%% + %.1f, %.1f%% + %.1f, %.2f + %.2f, %.2f + %.2f, %.0f%%, %.0f%%, %.2f + %.2f\n' ...
    , 'Coupled NMF', mean_err(17)*100, std_err(17)*100 ...
    , mean_err(18)*100, std_err(18)*100 ...
    , mean_err(19), std_err(19) ...
    , mean_err(20), std_err(20) ...      
    , mrna_consistency*100, protein_consistency*100 ...
    , mean_err(16), std_err(16));







