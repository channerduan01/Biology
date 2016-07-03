% There are 2 types of method:
% For 'MU' (default)
% Multiplicative update rules(Lee and Seung, 2001) and CNMF(Pauca et al., 2006)
%
% For 'ALS'
% Basic ALS(Paatero and Tapper, 1994)
function [W,H,iter,HIS] = mynmf(A,k,varargin)
cost = @(A,W,H) sqrt(sum(sum((A-W*H).^2)));

[m,n] = size(A);

par.method = 'MU';
par.alpha = 0;
par.beta = 0;
par.max_iter = 100;
par.min_iter = 20;
par.tol = 1e-3;
par.verbose = 0;

W = rand(m,k);
H = rand(k,n);
% only for test ===================================
for i = 1:n
    H(:,i) = H(:,i)/sum(H(:,i));
end
% =================================================

test_count = 0;

% Read optional parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'METHOD',              par.method = varargin{i+1};
            case 'ALPHA',               argAlpha = varargin{i+1};par.alpha = argAlpha;
            case 'BETA',                argBeta = varargin{i+1};par.beta = argBeta;
            case 'RATE',                par.rate = varargin{i+1};
            case 'W_INIT',              W = varargin{i+1};
            case 'H_INIT',              H = varargin{i+1};
            case 'MAX_ITER',            par.max_iter = varargin{i+1};
            case 'MIN_ITER',            par.min_iter = varargin{i+1};
            case 'TOL',                 par.tol = varargin{i+1};
            case 'VERBOSE',             par.verbose = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
if ~strcmp(par.method,'MU') && ~strcmp(par.method,'ALS') &&...
        ~strcmp(par.method,'ALS_W') && ~strcmp(par.method,'NMFSC')
    error('Unrecognized method: use ''MU'' or ''ALS'' or ''NMFSC''.');
end

if ~exist('argAlpha','var'), par.alpha = mean(A(:)); end
if ~exist('argBeta','var'), par.beta = mean(A(:)); end

par.start_cost = cost(A,W,H);
if par.verbose, display(par); end

initSC = getInitCriterion(A,W,H,par.alpha,par.beta);

% special operation for NMFSC
if strcmp(par.method,'NMFSC')
    for i = 1:length(H)
        H(:,i) = projection_operator(H(:,i),par.alpha,par.beta);
    end
end

SCconv = 0;
SC_COUNT = 3;
HIS = zeros(1,4);
for iter=1:par.max_iter
    switch par.method
        case 'MU'
            %     H = H .* (W'*A)./(W'*W*H + b*ones(r)*H + eps);
            H = H .* (W'*A)./(W'*W*H + par.beta*H + eps);
            % only for test ===================================
            for i = 1:n
                H(:,i) = H(:,i) + (1-sum(H(:,i)))/length(H(:,i));
            end
            H(H<0) = 0;
            % =================================================
            W = W .* (A*H')./(W*H*H' + par.alpha*W + eps);
        case 'ALS'
            H = pinv(W'*W+par.beta*eye(k))*(W'*A);
            H(H<0) = 0;
            W = (pinv(H*H'+par.alpha*eye(k))*(H*A'))';
            W(W<0) = 0;
        case 'ALS_W'
            H = pinv(W'*W+par.beta*eye(k))*(W'*A);
            W = (pinv(H*H'+par.alpha*eye(k))*(H*A'))';
            W(W<0) = 0;
        case 'NMFSC'
            % special step for sparsity matrix, a small rate is important
            %             H = H - par.rate*W'*(W*H-A); % every iteration just descend
            H = pinv(W'*W+0*eye(k))*(W'*A); % every iteration get best solution
            for i = 1:length(H)
                [H(:,i), ~] = projection_operator(H(:,i),par.alpha,par.beta);
                %                 test_count = test_count+iter
            end
            W = W .* (A*H')./(W*H*H' + eps);
    end
    SC = getStopCriterion(A,W,H,par.alpha,par.beta)/initSC;
    ver.iter = iter;
    ver.W_density = length(find(W>0))/(m*k);
    ver.H_density = length(find(H>0))/(n*k);
    ver.cost = cost(A,W,H);
    ver.SC = SC;
    HIS(iter,1) = ver.W_density;
    HIS(iter,2) = ver.H_density;
    HIS(iter,3) = ver.cost;
    HIS(iter,4) = ver.SC;
    if par.verbose, display(ver); end
    if (iter > par.min_iter)
        if SC <= par.tol
            SCconv = SCconv + 1;
            if (SCconv >= SC_COUNT), break; end
        else
            SCconv = 0;
        end
    end
end
end


%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------
function retVal = getStopCriterion(A,W,H,alpha,beta)
[gradW,gradH] = getGradient(A,W,H,alpha,beta);
pGrad = [gradW(gradW<0|W>0); gradH(gradH<0|H>0)];
pGradNorm = norm(pGrad);
retVal = pGradNorm/length(pGrad);
end
function retVal = getInitCriterion(A,W,H,alpha,beta)
[gradW,gradH] = getGradient(A,W,H,alpha,beta);
[m,k]=size(W);
[~,n]=size(H);
retVal = norm([gradW; gradH'],'fro')/((m*k)+(k*n));
end
function [gradW,gradH] = getGradient(A,W,H,alpha,beta)
gradW = W*(H*H') - A*H' + alpha*W;
gradH = (W'*W)*H - W'*A + beta*H;
end

