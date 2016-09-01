%%
% Tricky one, run all algorithms 

function [RESULT_TABLE, BEST_RESULT_TABLE, ERROR_NUM] ...
    = RunAllAlgorithm(REPEAT_NUM, MRNA, PROTEIN, THETA_ORIGINAL, H1_ORIGINAL, H2_ORIGINAL, K, J, T, N, bSHOW, ...
    thresholdRogers)
% basic init
max_iter = 100;
min_iter = 20;

[~, idx_mrna] = max(H1_ORIGINAL);
[~, idx_protein] = max(H2_ORIGINAL);

RESULT_TABLE = zeros(REPEAT_NUM, 4*5);
BEST_RESULT_TABLE = zeros(4, 3);
ERROR_NUM = zeros(4, 1);

IDX_MATRIX_MRNA = zeros(REPEAT_NUM, N, 4);
IDX_MATRIX_PROTEIN = zeros(REPEAT_NUM, N, 4);

% benchmark, double k-means
RESULT_KMEAN = cell(REPEAT_NUM,1);
for repeat_time = 1:REPEAT_NUM
    if bSHOW
    fprintf('\n k-means run on %d >>>>>>\n', repeat_time);
    end
    repeat_run = 1;
    if bSHOW
    opts = statset('Display','final');
    else
    opts = statset();
    end
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
    if bSHOW
    fprintf('k-means theta-error: %f, mrna correct rate: %f, protein correct rate: %f, total_cost: %f\n', ...
        theta_error, mrna_correct, protein_correct, total_cost);
    end
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

% benchmark, double NMF
RESULT_DNMF = cell(REPEAT_NUM,1);
repeat_time = 1;
error_result_num = 0;
while true
    if bSHOW
    fprintf('\n double NMF run on %d >>>>>>\n', repeat_time);
    end
    [W1,H1,W2,H2,~,HIS,last_iter] = ...
        CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
        , 'W_COEF', 0.2, 'H_COEF', 0.2 ...
        , 'VERBOSE', 0, 'METHOD', 'BP', 'PATIENCE', 0.1 ...
        , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
        );    
    total_cost = sum(HIS(last_iter+1,1:2));
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
        error_result_num = error_result_num + 1;
        if error_result_num > REPEAT_NUM
            throw(MException('Coupled NMF:failed','Too many errors!!!'));
        end
       continue; 
    end
    if bSHOW
    fprintf('Naive NMFs theta-error: %f, mrna correct rate: %f, protein correct rate: %f, total_cost: %f\n', ...
        theta_error, mrna_correct, protein_correct, total_cost);
    end
    RESULT_DNMF{repeat_time} = struct('total_cost',total_cost,'mrna_correct',mrna_correct, ...
        'protein_correct', protein_correct, 'theta_error',theta_error, 'theta', THETA2);
    RESULT_TABLE(repeat_time, 6:10) = [theta_error, mrna_correct, protein_correct, mean(sparsity(H1)), mean(sparsity(H2))];
    repeat_time = repeat_time + 1;
    if repeat_time > REPEAT_NUM, break; end
end
ERROR_NUM(2) = error_result_num;

best_dnmf_idx = 1;
for repeat_time = 2:REPEAT_NUM
    if RESULT_DNMF{repeat_time}.total_cost < RESULT_DNMF{best_dnmf_idx}.total_cost
        best_dnmf_idx = repeat_time;
    end
end
THETA2 = RESULT_DNMF{best_dnmf_idx}.theta;

% coupled Rogers
%----- for test
best_rogers_idx = 1;
%-----
[RESULT_, best_result, error_result_num] = ...
    MyCoupleClusteringRepeat(MRNA, PROTEIN, K, J, 100, thresholdRogers, bSHOW, REPEAT_NUM);
for repeat_time = 1:REPEAT_NUM
    [~, idx1] = max(RESULT_{repeat_time}.Q);
    [~, idx2] = max(RESULT_{repeat_time}.Q_J);
    IDX_MATRIX_MRNA(repeat_time,:,3) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,3) = idx2;
    THETA3 = ArrangeTheta(idx_mrna, idx1, idx_protein, idx2, RESULT_{repeat_time}.THETA, K, J);
    mrna_correct = purity(idx1, idx_mrna);
    protein_correct = purity(idx2, idx_protein);
    theta_error = norm(THETA_ORIGINAL-THETA3,'fro');
    if bSHOW
    fprintf('Rogers''s theta-error: %f, mrna correct rate: %f, protein correct rate: %f, low_bound: %f\n', ...
        theta_error, mrna_correct, protein_correct, RESULT_{repeat_time}.low_bound);
    fprintf('%f\n', mrna_correct+protein_correct);
    end
    RESULT_TABLE(repeat_time, 11:15) = [theta_error, mrna_correct, protein_correct, ...
        mean(sparsity(RESULT_{repeat_time}.Q)), mean(sparsity(RESULT_{repeat_time}.Q_J))];
end
ERROR_NUM(3) = error_result_num;

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


% coupled NMF
RESULT_CNMF = cell(REPEAT_NUM,1);
repeat_time = 1;
error_result_num = 0;
while true
    if bSHOW
    fprintf('\n Coupled NMF run on %d >>>>>>\n', repeat_time);
    end
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
        , 'W_COEF', 0.2, 'H_COEF', 0.2 ...
        , 'T_COEF', 1 ...
        , 'VERBOSE', 0, 'METHOD', 'BP', 'PATIENCE', 1 ...
        , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
        );
    [~, idx1] = max(H1);
    [~, idx2] = max(H2);
    IDX_MATRIX_MRNA(repeat_time,:,4) = idx1;
    IDX_MATRIX_PROTEIN(repeat_time,:,4) = idx2;
    total_cost = sum(HIS(last_iter+1,1:2));
    % recover THETA to original order!
    THETA4 = ArrangeTheta(idx_mrna_, idx1, idx_protein_, idx2, THETA4, K, J);
    mrna_correct = purity(idx1, idx_mrna_);
    protein_correct = purity(idx2, idx_protein_);
    theta_error = norm(THETA_ORIGINAL-THETA4,'fro');
    if isnan(total_cost) || sum(isnan(theta_error)) > 0
        error_result_num = error_result_num + 1;
        if error_result_num > REPEAT_NUM
            throw(MException('Coupled NMF:failed','Too many errors!!!'));
        end
       continue; 
    end
    if bSHOW
    fprintf('Coupled NMFs cost: %f, theta-error: %f, mrna correct rate: %f, protein correct rate: %f, total_cost: %f\n', ...
        cost_nmf, theta_error, mrna_correct, protein_correct, total_cost);
    end
    RESULT_CNMF{repeat_time} = struct('total_cost',total_cost,'mrna_correct',mrna_correct, ...
        'protein_correct', protein_correct, 'theta_error',theta_error, 'theta', THETA4);
    RESULT_TABLE(repeat_time, 16:20) = [theta_error, mrna_correct, protein_correct, mean(sparsity(H1)), mean(sparsity(H2))];
    repeat_time = repeat_time + 1;
    if repeat_time > REPEAT_NUM, break; end
end
ERROR_NUM(4) = error_result_num;

best_cnmf_idx = 1;
for repeat_time = 2:REPEAT_NUM
    if RESULT_CNMF{repeat_time}.total_cost < RESULT_CNMF{best_cnmf_idx}.total_cost
        best_cnmf_idx = repeat_time;
    end
end
THETA4 = RESULT_CNMF{best_cnmf_idx}.theta;


BEST_RESULT_TABLE(1, :) = [RESULT_TABLE(best_kmeans_idx, 2), RESULT_TABLE(best_kmeans_idx, 3), RESULT_TABLE(best_kmeans_idx, 1)];
BEST_RESULT_TABLE(2, :) = [RESULT_TABLE(best_dnmf_idx, 7), RESULT_TABLE(best_dnmf_idx, 8), RESULT_TABLE(best_dnmf_idx, 6)];
BEST_RESULT_TABLE(3, :) = [RESULT_TABLE(best_rogers_idx, 12), RESULT_TABLE(best_rogers_idx, 13), RESULT_TABLE(best_rogers_idx, 11)];
BEST_RESULT_TABLE(4, :) = [RESULT_TABLE(best_cnmf_idx, 17), RESULT_TABLE(best_cnmf_idx, 18), RESULT_TABLE(best_cnmf_idx, 16)];

end



