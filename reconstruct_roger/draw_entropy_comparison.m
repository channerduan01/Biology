%% init

close all
clear
clc

load RESULT
load RESULT_protein_permuted
load RESULT_mrna_protein_permuted

REPEAT = length(RESULT);


%% Result check
for i = 1:REPEAT
    if isnan(RESULT{i}.low_bound) 
        fprintf('nan in RESULT\n');
    end
    if isnan(RESULT_protein_permuted{i}.low_bound) 
        fprintf('nan in RESULT_protein_permuted\n');
    end
    
    if RESULT{i}.low_bound == -Inf || RESULT{i}.low_bound == +Inf
        fprintf('Inf in RESULT\n');
    end
    if RESULT_protein_permuted{i}.low_bound == -Inf || ...
            RESULT_protein_permuted{i}.low_bound == +Inf
        fprintf('Inf in RESULT_protein_permuted\n');
    end
end

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
    entropy_k_j_both_permuted(i) = RESULT_mrna_protein_permuted{i}.entropy_k_j;
    
    entropy_j_k_true(i) = RESULT{i}.entropy_j_k;
    entropy_j_k_permuted(i) = RESULT_protein_permuted{i}.entropy_j_k;
    entropy_j_k_both_permuted(i) = RESULT_mrna_protein_permuted{i}.entropy_j_k;
end

figure();
hold on
[pdft,x] = ksdensity(entropy_j_k_true);
plot(x,pdft)
[pdft,x] = ksdensity(entropy_j_k_permuted);
plot(x,pdft)
[pdft,x] = ksdensity(entropy_j_k_both_permuted);
plot(x,pdft)
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
[pdft,x] = ksdensity(entropy_k_j_both_permuted);
plot(x,pdft)
hold off
set(gca,'FontSize',20);
xlabel('Entropy of p(k | j)', 'FontSize', 16);
ylabel('Empirical density', 'FontSize', 16);
title('Entropy Analysis (100 runs)', 'FontSize', 20)
legend('true', 'protein permuted', 'both permuted');


%%



figure();
hold on
[pdft,x] = ksdensity(entropy_j_k_true);
plot(x,pdft)
[pdft,x] = ksdensity(entropy_j_k_permuted);
plot(x,pdft)
hold off
set(gca,'FontSize',20);
xlabel('Entropy of p(j | k)', 'FontSize', 16);
ylabel('Empirical density', 'FontSize', 16);
title('Entropy Analysis (100 runs)', 'FontSize', 20)
legend('true', 'protein permuted');









