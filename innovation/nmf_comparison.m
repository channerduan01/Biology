%% NMF clustering model
%
% init data
close all
% clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

% [MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] = GeneDataLoad();
% K = 15;
% J = 19;

REPEAT = 2;

IDX_MATRIX_MRNA = zeros(2, N);
IDX_MATRIX_PROTEIN = zeros(2, N);
K_test = K;
J_test = J;

index_ = 1;
[~, idx] = max(H1_ORIGINAL);
IDX_MATRIX_MRNA(index_,:) = idx;
[~, idx] = max(H2_ORIGINAL);
IDX_MATRIX_PROTEIN(index_,:) = idx;

for type_ = 1:2
    err_num = 0;
    i = 0;
    while true
        i = i + 1;
        if i > REPEAT, break; end
        try
            switch type_
                case 1
                    opts = statset('Display','final');
                    [idx_tmp, W1] = kmeans(MRNA',K,'Replicates',1,'Options',opts,'Start','sample');
                    W1 = W1';
                    H1 = zeros(K,length(idx_tmp));
                    for i = 1:length(idx_tmp)
                        H1(idx_tmp(i),i) = 1;
                    end
                    [idx_tmp, W2] = kmeans(MRNA',K,'Replicates',1,'Options',opts,'Start','sample');
                    W2 = W2';
                    H2 = zeros(K,length(idx_tmp));
                    for i = 1:length(idx_tmp)
                        H2(idx_tmp(i),i) = 1;
                    end             
                case 2
                    [W1,H1] = mynmf(MRNA,K_test,'METHOD','ALS','verbose',0,'MAX_ITER',100,'ALPHA',1,'BETA',1);
                    [W2,H2] = mynmf(PROTEIN,J_test,'METHOD','ALS','verbose',0,'MAX_ITER',100,'ALPHA',1,'BETA',1);
                case 3
                    [W1,H1] = mynmf(MRNA,K_test,'METHOD','ALS_W','verbose',0,'MAX_ITER',100,'ALPHA',1,'BETA',1);
                    [W2,H2] = mynmf(PROTEIN,J_test,'METHOD','ALS_W','verbose',0,'MAX_ITER',100,'ALPHA',1,'BETA',1);
                case 4
                    [W1,H1] = mynmf(MRNA,K_test,'verbose',0,'MAX_ITER',100,'ALPHA',1,'BETA',1);
                    [W2,H2] = mynmf(PROTEIN,J_test,'verbose',0,'MAX_ITER',100,'ALPHA',1,'BETA',1);          
                    %                     [W1,H1] = mynmf(MRNA,K_test,'METHOD','NMFSC','verbose',0,'MAX_ITER',10,'ALPHA',1,'BETA',1);
                    %                     [W2,H2] = mynmf(PROTEIN,J_test,'METHOD','NMFSC','verbose',0,'MAX_ITER',100,'ALPHA',1,'BETA',1);
                case 5
                    [W1,H1] = nmf(MRNA,K_test,'type','sparse','MAX_ITER',100,'verbose',0,'ALPHA',1,'BETA',1);
                    [W2,H2] = nmf(PROTEIN,J_test,'type','sparse','MAX_ITER',100,'verbose',0,'ALPHA',1,'BETA',1);
            end
        catch
            err_num = err_num+1;
            i = i - 1;
            continue;
        end
        
%         [~, idx] = max(H1);
%         IDX_MATRIX_MRNA(i,:) = idx;
%         [~, idx] = max(H2);
%         IDX_MATRIX_PROTEIN(i,:) = idx;
    end
    index_ = 2;
    [~, idx] = max(H1);
    IDX_MATRIX_MRNA(index_,:) = idx;
    [~, idx] = max(H2);
    IDX_MATRIX_PROTEIN(index_,:) = idx;
    [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
    [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
    fprintf('t-%d error: %d, mrna correct rate: %f, protein correct rate: %f\n', ...
        type_, err_num, mrna_consistency, protein_consistency);
    
%     [~, idx] = max(H1);
%     fprintf('For mRNA:\n');
%     for i = 1:K_test
%         fprintf('%d- %f%%\n', i, sum(idx==i)/N);
%     end
%     fprintf('For protein:\n');
%     [~, idx] = max(H2);
%     for i = 1:J_test
%         fprintf('%d- %f%%\n', i, sum(idx==i)/N);
%     end
end










