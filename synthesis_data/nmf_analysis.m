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


REPEAT_NUM = 5;
RESULT_TABLE = zeros(REPEAT_NUM, 4*3);
iter_used = zeros(REPEAT_NUM, 3);

IDX_MATRIX_MRNA = zeros(2, N);
IDX_MATRIX_PROTEIN = zeros(2, N);

index_ = 1;
[~, idx] = max(H1_ORIGINAL);
IDX_MATRIX_MRNA(index_,:) = idx;
[~, idx] = max(H2_ORIGINAL);
IDX_MATRIX_PROTEIN(index_,:) = idx;


for REPEAT_TIME = 1:REPEAT_NUM
fprintf('\n\n\nrun on %d >>>>>>\n', REPEAT_TIME);

%% bench mark, k-means
tic
index_ = 2;
repeat_run = 1;
opts = statset('Display','final');
[idx1,W1] = kmeans(MRNA',K,'Replicates',repeat_run,'Options',opts,'Start','sample');
W1 = W1';
IDX_MATRIX_MRNA(index_,:) = idx1;
H1 = zeros(K,length(idx1));
for i = 1:length(idx1)
    H1(idx1(i),i) = 1;
end
[idx2,W2] = kmeans(PROTEIN',J,'Replicates',repeat_run,'Options',opts,'Start','sample');
W2 = W2';
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
theta_error = norm(THETA_ORIGINAL-THETA1,'fro');
fprintf('k-means theta-error: %f, mrna correct rate: %f, protein correct rate: %f\n', ...
    theta_error, mrna_correct, protein_correct);
RESULT_TABLE(REPEAT_TIME, 1:3) = [theta_error, mrna_correct, protein_correct];
toc


%% SNMF innovative model, two stage
% while true
tic
index_ = 3;
[W1,H1,W2,H2,~,HIS,last_iter] = ...
    CoNMF_v2_separate(MRNA, PROTEIN, K, J ...
    , 'PATIENCE', 0.005 ...
    , 'W_COEF', w_coef, 'H_COEF', h_coef, 'T_COEF', t_coef ...
    , 'VERBOSE', 0, 'METHOD', method ...
    , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
    );
iter_used(REPEAT_TIME, 1) = iter_used(REPEAT_TIME, 1) + last_iter;
THETA2 = CalcuTheta(H1, H2, K, J, N);
[~, idx1] = max(H1);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(H2);
IDX_MATRIX_PROTEIN(index_,:) = idx2;



% recover THETA to original order!
THETA2 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA2, K, J);
[mrna_correct,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
theta_error = norm(THETA_ORIGINAL-THETA2,'fro');
fprintf('Naive NMFs theta-error: %f, mrna correct rate: %f, protein correct rate: %f\n', ...
    theta_error, mrna_correct, protein_correct);
RESULT_TABLE(REPEAT_TIME, 4:6) = [theta_error, mrna_correct, protein_correct];
toc
% if isnan(theta_error)
%    break; 
% end
% end

%% SNMF innovative model, coupled model
tic
index_ = 4;
[W1,H1,W2,H2,THETA3,HIS,last_iter] = ...
    CoNMF_v3_co2(MRNA, PROTEIN, K, J, 'MAX_ITER', 100, 'PATIENCE', 1 ...
    , 'W_COEF', w_coef, 'H_COEF', h_coef, 'T_COEF', t_coef ...
    , 'VERBOSE', 0, 'METHOD', method ...
    , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
    );
iter_used(REPEAT_TIME, 2) = iter_used(REPEAT_TIME, 2) + last_iter;
[~, idx1] = max(H1);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(H2);
IDX_MATRIX_PROTEIN(index_,:) = idx2;
% recover THETA to original order!
THETA3 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA3, K, J);
[mrna_correct,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
theta_error = norm(THETA_ORIGINAL-THETA3,'fro');
fprintf('Coupled NMFs theta-error: %f, mrna correct rate: %f, protein correct rate: %f\n', ...
    theta_error, mrna_correct, protein_correct);
RESULT_TABLE(REPEAT_TIME, 7:9) = [theta_error, mrna_correct, protein_correct];
toc

%% SNMF innovative model, flow model
tic
index_ = 5;
[W1,H1,W2,H2,THETA4,HIS,last_iter] = ...
    CoNMF_v4_flow(MRNA, PROTEIN, K, J ...
    , 'W_COEF', w_coef, 'H_COEF', h_coef, 'T_COEF', t_coef ...
    , 'VERBOSE', 0, 'METHOD', method ...
    , 'MIN_ITER', min_iter, 'MAX_ITER', max_iter ...
    );
iter_used(REPEAT_TIME, 3) = iter_used(REPEAT_TIME, 3) + last_iter;
[~, idx1] = max(H1);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(H2);
IDX_MATRIX_PROTEIN(index_,:) = idx2;
% recover THETA to original order!
THETA4 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA4, K, J);
[mrna_correct,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,index_],:));
[protein_correct,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,index_],:));
theta_error = norm(THETA_ORIGINAL-THETA4,'fro');
fprintf('FLow NMFs theta-error: %f, mrna correct rate: %f, protein correct rate: %f\n', ...
    theta_error, mrna_correct, protein_correct);
RESULT_TABLE(REPEAT_TIME, 10:12) = [theta_error, mrna_correct, protein_correct];
toc

end

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

%% visualize the result

% displayImages(W1, [20,20], 1, 'W1 result');






