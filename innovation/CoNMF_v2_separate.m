%%
% This is an intermediate product of the final algorithm
% This algorithm does not consider coupled connnection! So THETA does not
% work.
function [W1,H1,W2,H2,THETA,HIS,last_iter] = CoNMF_v2_separate(MRNA, PROTEIN, K, J, varargin)

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

% only for test, no use in this version of algorithm ===================================
THETA = rand(K,J);
for k = 1:K
    THETA(k,:) = THETA(k,:)/sum(THETA(k,:));
end
% ======================================================================================

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
record = calcuCost(V1,W1,H1,V2,W2,H2,THETA);
if par.verbose
    fprintf('CoNMF Init cost:%f  (V1-cost: %f, V2-cost: %f, Con-cost: %f - %f)\n', ...
        sum(record), record(1), record(2), record(3), record(4));
end
HIS(1, :) = record;
last_cost = sum(record);
for last_iter = 1:par.max_iter
    %     % plain mode, for test
    %     [H1,~,~] = nnlsm(W1,V1,H1);
    %     [W1,~,~] = nnlsm(H1',V1',W1');
    %     W1 = W1';
    %     [H2,~,~] = nnlsm(W2,V2,H2);
    %     [W2,~,~] = nnlsm(H2',V2',W2');
    %     W2 = W2';
    
    % regularized mode, real process
    [H1,~,~] = nnlsm([W1;eye_H_K],[V1;zeroKN],H1);
    [W1,~,~] = nnlsm([H1';eye_W_K],[V1';zeroKT],W1');
    W1 = W1';
    [H2,~,~] = nnlsm([W2;eye_H_J],[V2;zeroJN],H2);
    [W2,~,~] = nnlsm([H2';eye_W_J],[V2';zeroJT],W2');
    W2 = W2';
    
    record = calcuCost(V1,W1,H1,V2,W2,H2,THETA);
    new_cost = sum(record);
    HIS(last_iter+1, :) = record;
    step = new_cost - last_cost;
    last_cost = new_cost;
    if par.verbose
        fprintf('iter: %d, step: %f, cost: %f (V1-cost: %f, V2-cost: %f, Con-cost: %f & %f)\n', ...
            last_iter, step, new_cost, record(1), record(2), record(3), record(4));
    end
end
end
%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------
function [X,grad,iter] = nnlsm(A,B,init)
[X,grad,iter] = nnlsm_blockpivot(A,B,0,init);
end

function record = calcuCost(V1,W1,H1,V2,W2,H2,THETA)
cost = @(V,W,H) sqrt(sum(sum((V-W*H).^2)));
a = cost(V1,W1,H1);
b = cost(V2,W2,H2);
c = cost(W2,W1,THETA);
d = cost(W1,W2,pinv(THETA));
record = [a b c d];
end





