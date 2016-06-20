%% Roger EM coupled clustering model
%
% init data
close all
clear
clc

[MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] = GeneDataLoad();
K = 9;
J = 2;
MAX_ITER = 10;
REPEAT = 20;

%% taste
[W1,H1,W2,H2,THETA,HIS] = CoNMF(MRNA, PROTEIN, K, J, 5, 5, 3, MAX_ITER, true);
plot(HIS);
set(gca,'FontSize',16);
xlabel('iteration', 'FontSize', 16);
ylabel('cost', 'FontSize', 16);
title('Gradient', 'FontSize', 20)
legend('W1 and H1', 'W2 and H2', 'W1 and THETA', 'W2 and inverse THETA');

% %% Consistency
% IDX_MATRIX_MRNA = zeros(REPEAT, N);
% IDX_MATRIX_PROTEIN = zeros(REPEAT, N);
% for i = 1:REPEAT
%     [W1,H1,W2,H2,THETA,~] = CoNMF(MRNA, PROTEIN, K, J, 5, 5, 3, MAX_ITER, false);
%     [~, idx] = max(H1);
%     IDX_MATRIX_MRNA(i,:) = idx;
%     [~, idx] = max(H2);
%     IDX_MATRIX_PROTEIN(i,:) = idx;
% end
% [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
% [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
% fprintf('mrna consistency: %f, protein consistency: %f\n', mrna_consistency, protein_consistency);





