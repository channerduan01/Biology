%%
% This is an intermediate product of the final algorithm
% This algorithm does not consider coupled connnection! So THETA does not
% work.
function [W1,H1,W2,H2,THETA1,HIS,last_iter] = CoNMF_v3_co2(MRNA, PROTEIN, K, J, varargin)

% default values of parameters
par.method = '2STAGE';
par.wCoef = 1;
par.hCoef = 1;
par.tCoef = 1;
par.max_iter = 100;
par.min_iter = 10;
par.verbose = 0;

% init
[T,N] = size(MRNA);
V1 = MRNA;
V2 = PROTEIN;

W1 = rand(T,K);
H1 = rand(K,N);
W2 = rand(T,J);
H2 = rand(J,N);

% Key params, my coupled matrix
[C, THETA1, THETA2] = updateTheta(H1, H2, K, J, N);

% Read optional parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'METHOD',              par.method = varargin{i+1};
            case 'W_COEF',              par.wCoef = varargin{i+1};
            case 'H_COEF',              par.hCoef = varargin{i+1};
            case 'T_COEF',              par.tCoef = varargin{i+1};
            case 'W1_INIT',             W1 = varargin{i+1};
            case 'H1_INIT',             H1 = varargin{i+1};
            case 'W2_INIT',             W2 = varargin{i+1};
            case 'H2_INIT',             H2 = varargin{i+1};
            case 'MIN_ITER',            par.min_iter = varargin{i+1};
            case 'MAX_ITER',            par.max_iter = varargin{i+1};
            case 'VERBOSE',             par.verbose = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
if ~strcmp(par.method,'2STAGE') ...
        && ~strcmp(par.method,'4STAGE')
    error('Unrecognized method: use ''2STAGE'' or ''4STAGE''.');
end

% init basic parameters
eye_W_K = sqrt(par.wCoef)*eye(K);
eye_H_K = sqrt(par.hCoef)*eye(K);
eye_W_J = sqrt(par.wCoef)*eye(J);
eye_H_J = sqrt(par.hCoef)*eye(J);

zeroKN = zeros(K,N);
zeroKT = zeros(K,T);
zeroJN = zeros(J,N);
zeroJT = zeros(J,T);

% main process start
HIS = zeros(par.max_iter+1, 4);
record = calcuCost(V1,W1,H1,V2,W2,H2,THETA1,THETA2);
if par.verbose
    fprintf('CoNMF Init cost:%f  (V1-cost: %f, V2-cost: %f, Con-cost: %f - %f)\n', ...
        sum(record), record(1), record(2), record(3), record(4));
end
HIS(1, :) = record;
last_cost = sum(record);
for last_iter = 1:par.max_iter
    
    % first
    W1 = updateW(V1,W1,H1,par.wCoef);
    H2 = H2 + par.tCoef*(H1'*THETA1)';
%     H2 = H2.*(H1'*THETA1)';
    H2 = normalizeColumn(H2);
    W2 = updateW(V2,W2,H2,par.wCoef);
    
    % second
    H1 = updateH(V1,W1,H1,par.hCoef);
    H1 = normalizeColumn(H1);
    H2 = updateH(V2,W2,H2,par.hCoef);
    H2 = normalizeColumn(H2);
    [C, THETA1, THETA2] = updateTheta(H1, H2, K, J, N);
    
    % third
    W2 = updateW(V2,W2,H2,par.wCoef);
    H1 = H1 + par.tCoef*(H2'*THETA2)';
    W1 = updateW(V1,W1,H1,par.wCoef);
    
    % fourth
    H1 = updateH(V1,W1,H1,par.hCoef);
    H1 = normalizeColumn(H1);
    H2 = updateH(V2,W2,H2,par.hCoef);
    H2 = normalizeColumn(H2);
    [C, THETA1, THETA2] = updateTheta(H1, H2, K, J, N);    
    
    
    % check results
    record = calcuCost(V1,W1,H1,V2,W2,H2,THETA1,THETA2);
    new_cost = sum(record);
    HIS(last_iter+1, :) = record;
    step = last_cost - new_cost;
    last_cost = new_cost;
    if par.verbose
        fprintf('iter: %d, step: %f, cost: %f (V1-cost: %f, V2-cost: %f, Con-cost: %f & %f)\n', ...
            last_iter, step, new_cost, record(1), record(2), record(3), record(4));
    end
    if step > 0 && step < 1
        break;
    end
end
end
%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------
function [X,grad,iter] = nnlsm(A,B,init)
[X,grad,iter] = nnlsm_blockpivot(A,B,0,init);
end

function [C, tk, tj] = updateTheta(H1, H2, K, J, N)
C = H1*H2';
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
tk = tk ./ N;
tj = tj ./ N;
end

function record = calcuCost(V1,W1,H1,V2,W2,H2,THETA1,THETA2)
cost = @(V,W,H) sqrt(sum(sum((V-W*H).^2)));
a = cost(V1,W1,H1);
b = cost(V2,W2,H2);
c = norm((H1'*THETA1)'-H2,'fro');
d = norm((H2'*THETA2)'-H1,'fro');
record = [a b c d];
end

function A = normalizeRow(A)
for i = 1:size(A,1)
    A(i,:) = A(i,:)/sum(A(i,:));
end
A(A<0) = 0;
end

function A = normalizeColumn(A)
for i = 1:size(A,2)
    A(:,i) = A(:,i)/sum(A(:,i));
end
A(A<0) = 0;
end

% for test
function H = updateH(V,W,H,coef)
H = H .* (W'*V)./(W'*W*H + coef*H + eps);
end
function W = updateW(V,W,H,coef)
W = W .* (V*H')./(W*H*H' + coef*W + eps);
end

% function H = updateH(V,W,H,coef)
% H = pinv(W'*W+coef*eye(size(W,2)))*(W'*V);
% H(H<0) = 0;
% end
%
% function W = updateW(V,W,H,coef)
% W = (pinv(H*H'+coef*eye(size(H,1)))*(H*V'))';
% W(W<0) = 0;
% end


