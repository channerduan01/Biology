clc;
close all;
cost = @(V,W,H) sqrt(sum(sum((V-W*H).^2)));

% configuration
wCoef = 1;
hCoef = 2;
max_iter = 10; % should be small, because some algorithm may exit automatically

% data
V1 = MRNA;
W1 = rand(T,K);
H1 = rand(K,N);
V2 = PROTEIN;
W2 = rand(T,J);
H2 = rand(J,N);


TEST_NUM = 5;
HIS = zeros(TEST_NUM, 4);
for i = 1:TEST_NUM
    % Kim's result
    [W1_res,H1_res] = nmf(V1, K, 'type', 'regularized', 'nnls_solver', 'bp' ...
        , 'MAX_ITER',max_iter,'verbose', 1 ...
        , 'ALPHA', wCoef, 'BETA', hCoef ...
        ...%    , 'W_INIT', W1, 'H_INIT', H1 ...
        );
    [W2_res,H2_res] = nmf(V2, J, 'type', 'regularized', 'nnls_solver', 'bp' ...
        , 'MAX_ITER',max_iter,'verbose', 1 ...
        , 'ALPHA', wCoef, 'BETA', hCoef ...
        ...%    , 'W_INIT', W2, 'H_INIT', H2 ...
        );
    
    % My own result
    [W1_res1,H1_res1,W2_res1,H2_res1,THETA,~,~] = ...
        CoNMF_v2_separate(V1, V2, K, J ...
        , 'W_COEF', wCoef, 'H_COEF', hCoef, 'T_COEF', 1 ...
        , 'MAX_ITER', max_iter, 'verbose', 1 ...
        ...%    , 'W1_INIT', W1, 'H1_INIT', H1 ...
        ...%     , 'W2_INIT', W2, 'H2_INIT', H2 ...
        );
    
    HIS(i, :) = [cost(V1,W1_res,H1_res), cost(V2,W2_res,H2_res), ...
        cost(V1,W1_res1,H1_res1), cost(V2,W2_res1,H2_res1)];
end

% Output
for i = 1:TEST_NUM
    fprintf('Kim''s tool: %f %f\nMy results: %f %f\n', HIS(i, 1), HIS(i, 2), HIS(i, 3), HIS(i, 4));
end
%% draw the distribution
figure();
hold on
[pdft,x] = ksdensity(HIS(:,1));
plot(x,pdft)
[pdft,x] = ksdensity(HIS(:,2));
plot(x,pdft)
[pdft,x] = ksdensity(HIS(:,3));
plot(x,pdft)
[pdft,x] = ksdensity(HIS(:,4));
plot(x,pdft)
hold off
set(gca,'FontSize',20);
xlabel('cost', 'FontSize', 16);
ylabel('Empirical density', 'FontSize', 16);
title('Algorithm verification', 'FontSize', 20)
legend('Kim''s V1', 'Kim''s V2', 'My V1', 'My V2');


