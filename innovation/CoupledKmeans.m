%% This algorithm does not work~
% 
%
function [W1,H1,W2,H2,THETA1,cost_] = CoupledKmeans(V1, V2, K, J, interact_rate, MAX_ITER)
cost = @(V,W,H) norm(V-W*H,'fro');
N = size(V1, 2);
% init
ii = randperm(N);
W1 = V1(:,ii(1:K));
H1 = zeros(K,N);
ii = randperm(N);
W2 = V2(:,ii(1:J));
H2 = zeros(J,N);
cost_ = cost(V1,W1,H1)+cost(V2,W2,H2);
fprintf('initial cost: %f\n', cost_);
for iter = 1:MAX_ITER
    % expected
    H1 = zeros(K,N);
    H2 = zeros(J,N);
    dists = zeros(1,N);
    for i = 1:N
        [dists(i),idx_min] = min(sum((W1-repmat(V1(:,i),1,K)).^2));
        H1(idx_min,i) = 1;
        [dists(i),idx_min] = min(sum((W2-repmat(V2(:,i),1,J)).^2));
        H2(idx_min,i) = 1;        
    end
    % interact
    [~, THETA1, THETA2] = calcuCorrelation(H1,H2,K,J,N);
    WaveForH2 = (H1'*THETA1)';
    WaveForH1 = (H2'*THETA2)';
    WaveForH2 = WaveForH2.*H2;
    WaveForH1 = WaveForH1.*H1;
    H2_ = (1-interact_rate)*H2 + interact_rate*WaveForH2;
    H2 = zeros(size(H2_));
    H2(H2_>0) = 1; 
    H1_ = (1-interact_rate)*H1 + interact_rate*WaveForH1;
    H1 = zeros(size(H1_));
    H1(H1_>0) = 1;
    % maximize
    for i = 1:K
        W1(:,i) = mean(V1(:,H1(i,:)==1),2);
    end
    for i = 1:J
        W2(:,i) = mean(V2(:,H2(i,:)==1),2);
    end
    if sum(any(isnan(W1))) > 0 || sum(any(isnan(W2))) > 0
        error('MATLAB:nan cluster centre!');
    end
    cost_ = cost(V1,W1,H1)+cost(V2,W2,H2);
    fprintf('%d cost: %f\n', iter, cost_);
end
fprintf('total cost: %f, total-dists: %f\n', cost_, sum(dists));
end



function [C, tk, tj] = calcuCorrelation(H1, H2, K, J, N)
C = H1*H2';
C = C ./ N;
tk = zeros(K,J);
tj = zeros(J,K);
PI_K_ = sum(H1,2)/N;
PI_J_ = sum(H2,2)/N;
for j = 1:J
    tk(:,j) = C(:,j)./PI_K_;
end
for k = 1:K
    tj(:,k) = C(k,:)./PI_J_';
end
tk(isnan(tk)) = 0;
tj(isnan(tj)) = 0;
end