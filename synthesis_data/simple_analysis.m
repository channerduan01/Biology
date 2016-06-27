%% 
% Methods comparison

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

IDX_MATRIX_MRNA = zeros(2, N);
IDX_MATRIX_PROTEIN = zeros(2, N);

index_ = 1;
[~, idx] = max(H1_ORIGINAL);
IDX_MATRIX_MRNA(index_,:) = idx;
[~, idx] = max(H2_ORIGINAL);
IDX_MATRIX_PROTEIN(index_,:) = idx;


%% bench mark, k-means
index_ = 2;
repeat = 1;
opts = statset('Display','final');
[idx1,~] = kmeans(MRNA',K,'Replicates',repeat,'Options',opts,'Start','sample');
IDX_MATRIX_MRNA(index_,:) = idx1;
H1 = zeros(K,length(idx1));
for i = 1:length(idx1)
    H1(idx1(i),i) = 1;
end
[idx2,~] = kmeans(PROTEIN',J,'Replicates',repeat,'Options',opts,'Start','sample');
IDX_MATRIX_PROTEIN(index_,:) = idx2;
H2 = zeros(J,length(idx2));
for i = 1:length(idx2)
    H2(idx2(i),i) = 1;
end
THETA1 = CalcuTheta(H1, H2, K, J, N);
% recover THETA to original order!
THETA1 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA1, K, J);

[mrna_consistency,mrna_C1] = CalcuConsistency(IDX_MATRIX_MRNA([1,2],:));
[protein_consistency,protein_C1] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,2],:));
fprintf('k-means mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);
fprintf('k-means norm-fro of gap: %.2f, norm-1 of gap: %.2f\n', ...
    norm(THETA1-THETA_ORIGINAL,'fro'),norm(THETA1-THETA_ORIGINAL,1));

% %% SNMF innovative model
% index_ = 3;
% [W1,H1,W2,H2,THETA2,HIS] = ...
%     CoNMF(MRNA, PROTEIN, K, J, 0.1, 0.01, 0.01, 100, true, 0.1);
% [~, idx] = max(H1);
% IDX_MATRIX_MRNA(index_,:) = idx;
% [~, idx] = max(H2);
% IDX_MATRIX_PROTEIN(index_,:) = idx;
% for k = 1:K
%     THETA2(k,:) = THETA2(k,:)./sum(THETA2(k,:));
% end
% 
% [mrna_consistency,mrna_C2] = CalcuConsistency(IDX_MATRIX_MRNA([1,3],:));
% [protein_consistency,protein_C2] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,3],:));
% fprintf('SNMF mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);
% fprintf('SNMF norm-fro of gap: %.2f, norm-1 of gap: %.2f\n', ...
%     norm(THETA2-THETA_ORIGINAL,'fro'),norm(THETA2-THETA_ORIGINAL,1));

%% Roger's statistic model
index_ = 4;
[Q,R,PI_K,AVG_K,VARIANCE_K,THETA3,AVG_J,VARIANCE_J] = ...
    MyCoupleClustering(MRNA, PROTEIN, K, J, 100, 1, true);
R_J = zeros(J, N);
for j = 1:J
    for i = 1:N
        tmp = R(j,i,:);
        R_J(j,i) = sum(tmp(:).*PI_K);
    end
end
[~, idx1] = max(Q);
IDX_MATRIX_MRNA(index_,:) = idx1;
[~, idx2] = max(R_J);
IDX_MATRIX_PROTEIN(index_,:) = idx2;

% recover THETA to original order!
THETA3 = ArrangeTheta(IDX_MATRIX_MRNA(1,:), idx1, IDX_MATRIX_PROTEIN(1,:), idx2, THETA3, K, J);


[mrna_consistency,mrna_C3] = CalcuConsistency(IDX_MATRIX_MRNA([1,4],:));
[protein_consistency,protein_C3] = CalcuConsistency(IDX_MATRIX_PROTEIN([1,4],:));
fprintf('Roger''s mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);
fprintf('Roger''s norm-fro of gap: %.2f, norm-1 of gap: %.2f\n', ...
    norm(THETA3-THETA_ORIGINAL,'fro'),norm(THETA3-THETA_ORIGINAL,1));


%% compare the patterns
figure();
hold on;
subplot(221), imagesc(THETA_ORIGINAL);
subplot(222), imagesc(THETA1);
% subplot(223), imagesc(THETA2);
subplot(224), imagesc(THETA3);
hold off;



