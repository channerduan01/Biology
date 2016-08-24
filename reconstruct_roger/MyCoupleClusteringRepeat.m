function [RESULT, best_result, error_result_num] = ...
    MyCoupleClusteringRepeat(MRNA, PROTEIN, K, J, MAX_ITER, patience, b_verbose, REPEAT)
[T1, N] = size(MRNA);
[T2, ~] = size(PROTEIN);
RESULT = cell(REPEAT,1);
index = 0;
error_result_num = 0;
while true
    % -------- permuted Proteins, keep annotation
    %     ii = randperm(N);
    %     MRNA(:,1:N) = MRNA(:,ii);
    %     PROTEIN(:,1:N) = PROTEIN(:,ii);
    % --------
    index = index + 1;
    if index > REPEAT, break;end
    fprintf('\nstart iter: %d >>>>>>\n', index);
    try
        [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = ...
            MyCoupleClustering(MRNA, PROTEIN, K, J, MAX_ITER, patience, b_verbose);
        [THETA_reverse, entropy_j_k, entropy_k_j] = EntropyCalculate(K, J, PI_K, THETA);
        [R_J, Q_J] = CalcuSubclusterBelonging(MRNA, AVG_K, VARIANCE_K, PROTEIN, AVG_J, VARIANCE_J, PI_K, THETA_reverse, R, K, J, N, T1, T2);
        low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T1,T2,N);
    catch ME
        low_bound = NaN;
        fprintf('error: %s\n', getReport(ME))
    end
    if isnan(low_bound) || low_bound == -Inf || low_bound == Inf
        error_result_num = error_result_num + 1;
        if error_result_num > 10*REPEAT
            throw(MException('MyCoupleClusteringRepeat:failed','Too many errors!!!'));
        end
        index = index - 1;
        continue;
    end
    RESULT{index} = struct('low_bound',low_bound,'entropy_j_k',entropy_j_k,'entropy_k_j',entropy_k_j, ...
        'THETA_reverse',THETA_reverse,'Q',Q,'R',R,'R_J',R_J,'Q_J',Q_J,'PI_K',PI_K,'AVG_K',AVG_K, ...
        'VARIANCE_K',VARIANCE_K,'THETA',THETA,'AVG_J',AVG_J,'VARIANCE_J',VARIANCE_J);
end
best_result = RESULT{1};
for i = 2:REPEAT
    if best_result.low_bound < RESULT{i}.low_bound
        best_result = RESULT{i};
    end
end




