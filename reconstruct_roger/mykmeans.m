function [W,H] = mykmeans(V, K, MAX_ITER)
cost = @(V,W,H) norm(V-W*H,'fro');

N = size(V, 2);
% init
ii = randperm(N);
W = V(:,ii(1:K));
H = zeros(K,N);
for iter = 1:MAX_ITER
    cost_ = cost(V,W,H);
    fprintf('%d cost: %f\n', iter, cost_);
    % expected
    H = zeros(K,N);
    dists = zeros(1,N);
    for i = 1:N
        [dists(i),idx_min] = min(sum((W-repmat(V(:,i),1,K)).^2));
        H(idx_min,i) = 1;
    end
    
    % maximize
    for i = 1:K
        W(:,i) = mean(V(:,H(i,:)==1),2);
    end
    if sum(any(isnan(W))) > 0
        error('MATLAB:nan cluster centre!');
    end
end
fprintf('total cost: %f, total-dists: %f\n', cost_, sum(dists));

