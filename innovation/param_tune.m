%% NMF coupled clustering model
%
% init data
close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

load hmec;
MRNA = data{1};
P1 = data{2};
MRNA = MRNA';
PROTEIN = P1(:,2:7)';
PROTEIN_ORIGINAL = P1';
MRNA = normalize(MRNA);
PROTEIN = normalize(PROTEIN);
[T,N] = size(MRNA);

%s basic params
REPEAT = 20;
K = 15;
J = 19;
min_iter = 20;
max_iter = 200;
switch_on_verbose = false;
    
    
% param search field
W_coef = [0.1];
H_coef = [0.1];
T_coef = [0.5];
pattern_num = length(W_coef)*length(H_coef)*length(T_coef);
RESULT_PARAM = cell(pattern_num,1);
index_ = 1;
for i1 = 1:length(W_coef)
    for i2 = 1:length(H_coef)
        for i3 = 1:length(T_coef)
            RESULT_PARAM{index_} = struct('wCoef', W_coef(i1),'hCoef', H_coef(i2), ...
                'tCoef', T_coef(i3), 'pattern', strcat(num2str(W_coef(i1)),'-',num2str(H_coef(i2)),'-',num2str(T_coef(i2)) ));
        end
    end
end

for pa = 1:pattern_num
    fprintf('\n\n >>>>>> pattern: %d <<<<<<< %s\n', pa, RESULT_PARAM{pa}.pattern);
    RESULT = cell(REPEAT,1);
    IDX_MATRIX_MRNA = zeros(REPEAT, N);
    IDX_MATRIX_PROTEIN = zeros(REPEAT, N);
    err_num = 0;
    i = 0;
    HIS = zeros(REPEAT, 4);
    
    while true
        i = i + 1;
        if i > REPEAT, break; end
        fprintf('process-iter >> %d\n', i);
        [W1_res,H1_res,W2_res,H2_res,Theta,HIS,last_iter] = CoNMF_v4_flow(MRNA, PROTEIN, K, J ...
            , 'MAX_ITER', max_iter, 'MIN_ITER', min_iter, 'VERBOSE', switch_on_verbose, 'METHOD', 'BP' ...
            , 'W_COEF', RESULT_PARAM{pa}.wCoef, 'H_COEF', RESULT_PARAM{pa}.hCoef, 'T_COEF', RESULT_PARAM{pa}.tCoef, 'PATIENCE', 0.0001 ...
            );
        RESULT{i} = struct('W1',W1_res,'H1',H1_res,'W2',W2_res,'H2',H2_res, ...
            'Theta',Theta,'HIS',HIS,'last_iter',last_iter);
        
        if sum(sum(isnan(Theta))) > 0
            err_num = err_num+1;
            i = i - 1;
            continue;
        end
        [~, idx] = max(H1_res);
        IDX_MATRIX_MRNA(i,:) = idx;
        [~, idx] = max(H2_res);
        IDX_MATRIX_PROTEIN(i,:) = idx;
    end
    
    if REPEAT > 1
        [mrna_consistency,mrna_C] = CalcuConsistency(IDX_MATRIX_MRNA);
        [protein_consistency,protein_C] = CalcuConsistency(IDX_MATRIX_PROTEIN);
        fprintf('error: %d, mrna consistency: %f, protein consistency: %f\n', ...
            err_num, mrna_consistency, protein_consistency);
    else
        mrna_consistency = 0;
        protein_consistency = 0;
    end
    RESULT_PARAM{pa}.err_num = err_num;
    RESULT_PARAM{pa}.mrna_con = mrna_consistency;
    RESULT_PARAM{pa}.protein_con = protein_consistency;
    
    ALL_HIS = zeros(max_iter+1,4);
    iters = zeros(1,REPEAT);
    for i = 1:REPEAT
        ALL_HIS = ALL_HIS + RESULT{i}.HIS;
        iters(1,i) = RESULT{i}.last_iter;
    end
    ALL_HIS = ALL_HIS/REPEAT;
    RESULT_PARAM{pa}.HIS = ALL_HIS;
    RESULT_PARAM{pa}.ITER = [num2str(mean(iters)),'+',num2str(std(iters))];
end

%%
















