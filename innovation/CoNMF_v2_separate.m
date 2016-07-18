%%
% This is an intermediate product of the final algorithm
% This algorithm does not consider coupled connnection! So THETA does not
% work.
function [W1,H1,W2,H2,THETA,HIS,last_iter] = CoNMF_v2_separate(MRNA, PROTEIN, K, J, varargin)

% default values of parameters
par.method = 'MU';
par.wCoef = 1;
par.hCoef = 1;
par.tCoef = 1;
par.max_iter = 100;
par.min_iter = 10;
par.verbose = 0;
par.patience = 0.1;

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
            case 'PATIENCE',            par.patience = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
if ~strcmp(par.method,'MU') ...
        && ~strcmp(par.method,'ALS') ...
        && ~strcmp(par.method,'BP') ...
        && ~strcmp(par.method,'SBP')
    error('Unrecognized method: use ''MU'' or ''ALS'' or ''BP'' or ''SBP''.');
end

% init basic parameters
par.eye_W_K = sqrt(par.wCoef)*eye(K);
par.eye_H_K = sqrt(par.hCoef)*eye(K);
par.eye_W_J = sqrt(par.wCoef)*eye(J);
par.eye_H_J = sqrt(par.hCoef)*eye(J);

par.zeroKN = zeros(K,N);
par.zeroKT = zeros(K,T);
par.zeroJN = zeros(J,N);
par.zeroJT = zeros(J,T);
par.zero1N = zeros(1,N);

par.ones_H_K = sqrt(par.hCoef)*ones(1,K);
par.ones_H_J = sqrt(par.hCoef)*ones(1,J);


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
    
    
    H1 = updateH(V1,W1,H1,par,1);
    H2 = updateH(V2,W2,H2,par,2);
    
    W1 = updateW(V1,W1,H1,par,1);
    W2 = updateW(V2,W2,H2,par,2);
    
    
    record = calcuCost(V1,W1,H1,V2,W2,H2,THETA);
    new_cost = sum(record);
    HIS(last_iter+1, :) = record;
    step = last_cost - new_cost;
    last_cost = new_cost;
    if par.verbose
        fprintf('iter: %d, step: %f, cost: %f (V1-cost: %f, V2-cost: %f, Con-cost: %f & %f)\n', ...
            last_iter, step, new_cost, record(1), record(2), record(3), record(4));
    end
    if step >= 0 && step < par.patience
        break;
    end
end
end
%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------

function record = calcuCost(V1,W1,H1,V2,W2,H2,THETA)
cost = @(V,W,H) sqrt(sum(sum((V-W*H).^2)));
a = cost(V1,W1,H1);
b = cost(V2,W2,H2);
c = 0;
d = 0;
record = [a b c d];
end

function H = updateH(V,W,H,par,idx)
if strcmp(par.method,'MU')
    H = updateH_MU(V,W,H,par.hCoef);
elseif strcmp(par.method,'ALS')
    H = updateH_ALS(V,W,H,par.hCoef);
elseif strcmp(par.method,'BP')
    if idx == 1
        [H,~,~] = nnlsm([W;par.eye_H_K],[V;par.zeroKN],H);
    elseif idx == 2
        [H,~,~] = nnlsm([W;par.eye_H_J],[V;par.zeroJN],H);
    else
        error(['updateH param idx unavailable: ', idx]);
    end
elseif strcmp(par.method,'SBP')
    if idx == 1
        [H,~,~] = nnlsm([W;par.ones_H_K],[V;par.zero1N],H);
    elseif idx == 2
        [H,~,~] = nnlsm([W;par.ones_H_J],[V;par.zero1N],H);
    else
        error(['updateH param idx unavailable: ', idx]);
    end
end
end

function W = updateW(V,W,H,par,idx)
if strcmp(par.method,'MU')
    W = updateW_MU(V,W,H,par.hCoef);
elseif strcmp(par.method,'ALS')
    W = updateW_ALS(V,W,H,par.hCoef);
elseif strcmp(par.method,'BP')
    if idx == 1
        [W,~,~] = nnlsm([H';par.eye_W_K],[V';par.zeroKT],W');
        W = W';
    elseif idx == 2
        [W,~,~] = nnlsm([H';par.eye_W_J],[V';par.zeroJT],W');
        W = W';
    else
        error(['updateW param idx unavailable: ', idx]);
    end
elseif strcmp(par.method,'SBP') % absolutely
    if idx == 1
        [W,~,~] = nnlsm([H';par.eye_W_K],[V';par.zeroKT],W');
        W = W';
    elseif idx == 2
        [W,~,~] = nnlsm([H';par.eye_W_J],[V';par.zeroJT],W');
        W = W';
    else
        error(['updateW param idx unavailable: ', idx]);
    end
end
end

function [X,grad,iter] = nnlsm(A,B,init)
[X,grad,iter] = nnlsm_blockpivot(A,B,0,init);
end

function H = updateH_MU(V,W,H,coef)
H = H .* (W'*V)./(W'*W*H + coef*H + eps);
end
function W = updateW_MU(V,W,H,coef)
W = W .* (V*H')./(W*H*H' + coef*W + eps);
end

function H = updateH_ALS(V,W,H,coef)
H = pinv(W'*W+coef*eye(size(W,2)))*(W'*V);
H(H<0) = 0;
end

function W = updateW_ALS(V,W,H,coef)
W = (pinv(H*H'+coef*eye(size(H,1)))*(H*V'))';
W(W<0) = 0;
end


