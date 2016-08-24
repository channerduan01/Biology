%% init

close all
clear
clc

% load RESULT
% load RESULT_protein_permuted
% load RESULT_mrna_protein_permuted

load RESULT_
load RESULT_protein_permuted_

% load NMF_RESULT

REPEAT = length(RESULT);

% %% Result check
% for i = 1:REPEAT
%     if isnan(RESULT{i}.low_bound) 
%         fprintf('nan in RESULT\n');
%     end
%     if isnan(RESULT_protein_permuted{i}.low_bound) 
%         fprintf('nan in RESULT_protein_permuted\n');
%     end
%     
%     if RESULT{i}.low_bound == -Inf || RESULT{i}.low_bound == +Inf
%         fprintf('Inf in RESULT\n');
%     end
%     if RESULT_protein_permuted{i}.low_bound == -Inf || ...
%             RESULT_protein_permuted{i}.low_bound == +Inf
%         fprintf('Inf in RESULT_protein_permuted\n');
%     end
% end

%% Sparsity calculate for NMF
res_ = zeros(REPEAT, 4);
for i = 1:REPEAT
    res_(i,1) = mean(sparsity(RESULT{i}.H1));
    res_(i,2) = mean(sparsity(RESULT{i}.H2));
    res_(i,3) = mean(sparsity(RESULT_protein_permuted{i}.H1));
    res_(i,4) = mean(sparsity(RESULT_protein_permuted{i}.H2));    
end

%% Sparsity calculate for Rogers
res_ = zeros(REPEAT, 4);
for i = 1:REPEAT
    res_(i,1) = mean(sparsity(RESULT{i}.Q));
    res_(i,2) = mean(sparsity(RESULT{i}.Q_J));
    res_(i,3) = mean(sparsity(RESULT_protein_permuted{i}.Q));
    res_(i,4) = mean(sparsity(RESULT_protein_permuted{i}.Q_J));    
end

%% Consistency calculate
% N = length(RESULT{1}.H1);
N = length(RESULT{1}.Q);
IDX_MATRIX_MRNA = zeros(REPEAT, N);
IDX_MATRIX_PROTEIN = zeros(REPEAT, N);
for i = 1:REPEAT
%     [~, idx] = max(RESULT{i}.H1);
    [~, idx] = max(RESULT{i}.Q);
    IDX_MATRIX_MRNA(i,:) = idx;
%     [~, idx] = max(RESULT{i}.H2);
    [~, idx] = max(RESULT{i}.Q_J);
    IDX_MATRIX_PROTEIN(i,:) = idx;    
end
[mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
[protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);

%% empirical
entropy_k_j_true = zeros(1,REPEAT);
entropy_k_j_permuted = zeros(1,REPEAT);
entropy_k_j_both_permuted = zeros(1,REPEAT);
entropy_j_k_true = zeros(1,REPEAT);
entropy_j_k_permuted = zeros(1,REPEAT);
entropy_j_k_both_permuted = zeros(1,REPEAT);
for i = 1:REPEAT
    entropy_k_j_true(i) = RESULT{i}.entropy_k_j;
    entropy_k_j_permuted(i) = RESULT_protein_permuted{i}.entropy_k_j;
%     entropy_k_j_both_permuted(i) = RESULT_mrna_protein_permuted{i}.entropy_k_j;
    
    entropy_j_k_true(i) = RESULT{i}.entropy_j_k;
    entropy_j_k_permuted(i) = RESULT_protein_permuted{i}.entropy_j_k;
%     entropy_j_k_both_permuted(i) = RESULT_mrna_protein_permuted{i}.entropy_j_k;
end

figure();
hold on
[pdft,x] = ksdensity(entropy_j_k_true);
plot(x,pdft)
[pdft,x] = ksdensity(entropy_j_k_permuted);
plot(x,pdft)
% [pdft,x] = ksdensity(entropy_j_k_both_permuted);
% plot(x,pdft)
hold off
set(gca,'FontSize',20);
xlabel('Entropy of p(j | k)', 'FontSize', 16);
ylabel('Empirical density', 'FontSize', 16);
title('Entropy Analysis (100 runs)', 'FontSize', 20)
legend('true', 'protein permuted', 'both permuted');


figure();
hold on
[pdft,x] = ksdensity(entropy_k_j_true);
plot(x,pdft)
[pdft,x] = ksdensity(entropy_k_j_permuted);
plot(x,pdft)
% [pdft,x] = ksdensity(entropy_k_j_both_permuted);
% plot(x,pdft)
hold off
set(gca,'FontSize',20);
xlabel('Entropy of p(k | j)', 'FontSize', 16);
ylabel('Empirical density', 'FontSize', 16);
title('Entropy Analysis (100 runs)', 'FontSize', 20)
legend('true', 'protein permuted', 'both permuted');


%%



% figure();
% hold on
% [pdft,x] = ksdensity(entropy_j_k_true);
% plot(x,pdft)
% [pdft,x] = ksdensity(entropy_j_k_permuted);
% plot(x,pdft)
% hold off
% set(gca,'FontSize',20);
% xlabel('Entropy of p(j | k)', 'FontSize', 16);
% ylabel('Empirical density', 'FontSize', 16);
% title('Entropy Analysis (100 runs)', 'FontSize', 20)
% legend('true', 'protein permuted');









