addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

IDX_MATRIX_MRNA = zeros(2, N);
IDX_MATRIX_PROTEIN = zeros(2, N);

i = 1;
[~, idx] = max(H1_ORIGINAL);
IDX_MATRIX_MRNA(i,:) = idx;
[~, idx] = max(H2_ORIGINAL);
IDX_MATRIX_PROTEIN(i,:) = idx;


i = 2;
[W1,H1,W2,H2,THETA,HIS] = ...
    CoNMF(MRNA, PROTEIN, K, J, 0.1, 0.01, 0.01, 100, true, 0.1);
[~, idx] = max(H1);
IDX_MATRIX_MRNA(i,:) = idx;
[~, idx] = max(H2);
IDX_MATRIX_PROTEIN(i,:) = idx;

% [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = ...
%     MyCoupleClustering(MRNA, PROTEIN, K, J, 100, 50, true);
% R_J = zeros(J, N);
% for j = 1:J
%     for i = 1:N
%         tmp = R(j,i,:);
%         R_J(j,i) = sum(tmp(:).*PI_K);
%     end
% end
% [~, idx] = max(Q);
% IDX_MATRIX_MRNA(i,:) = idx;
% [~, idx] = max(R_J);
% IDX_MATRIX_PROTEIN(i,:) = idx;



%%

[mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
[protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
fprintf('mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);





