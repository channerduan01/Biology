%%
% Data Generator
%

function [MRNA_ORIGINAL, PROTEIN_ORIGINAL, H1_ORIGINAL, H2_ORIGINAL, AVG_K_, COV_K_, AVG_J_, COV_J_] ...
    = GenerateBasisData(K, J, T, N)
VALUE_RANGE_D = 1:0.1:2;
VALUE_VARIANCE = 0.3;

OVERLAP_COEFF = 1;

% Init distributions
while true
    AVG_K_ = zeros(T, K);
    COV_K_ = zeros(T, T, K);
    for extra = 1:OVERLAP_COEFF
        for i = 1:T
            k = floor(rand()*K)+1;
            AVG_K_(i,k) = VALUE_RANGE_D(floor(rand()*length(VALUE_RANGE_D))+1);
            COV_K_(i,i,k) = VALUE_VARIANCE;
        end
    end
    pass_ = true;
    for k = 1:K
        if sum(AVG_K_(:,k)) == 0
           pass_ = false;
           break; 
        end
    end
    if pass_ 
        break;
    end
end
while true
    AVG_J_ = zeros(T, J);
    COV_J_ = zeros(T, T, J);    
    for extra = 1:OVERLAP_COEFF
        for i = 1:T
            j = floor(rand()*J)+1;
            AVG_J_(i,j) = VALUE_RANGE_D(floor(rand()*length(VALUE_RANGE_D))+1);
            COV_J_(i,i,j) = VALUE_VARIANCE;
        end
    end
    pass_ = true;
    for j = 1:J
        if sum(AVG_J_(:,j)) == 0
           pass_ = false;
           break; 
        end
    end
    if pass_ 
        break;
    end
end
% generate original data
MRNA_ORIGINAL = zeros(T, N);
PROTEIN_ORIGINAL = zeros(T, N);
H1_ORIGINAL = zeros(K, N);
H2_ORIGINAL = zeros(J, N);
range = @(K,N,k) floor(N/K*(k-1))+1:floor(N/K*k);
J_RANGE = randperm(K);
for k = 1:K
    k_ = J_RANGE(k);
    H1_ORIGINAL(k, range(K,N,k)) = 1;
    H2_ORIGINAL(k_, range(K,N,k)) = 1;
    MRNA_ORIGINAL(:, range(K,N,k)) = mvnrnd(AVG_K_(:,k)', COV_K_(:,:,k), length(range(K,N,k)))';
    PROTEIN_ORIGINAL(:, range(K,N,k)) = mvnrnd(AVG_J_(:,k_)', COV_J_(:,:,k_), length(range(K,N,k)))';    
end
ii = randperm(N);
MRNA_ORIGINAL = MRNA_ORIGINAL(:,ii);
PROTEIN_ORIGINAL = PROTEIN_ORIGINAL(:,ii);
H1_ORIGINAL = H1_ORIGINAL(:,ii);
H2_ORIGINAL = H2_ORIGINAL(:,ii);

% MRNA_ORIGINAL(MRNA_ORIGINAL<0) = 0;
% PROTEIN_ORIGINAL(PROTEIN_ORIGINAL<0) = 0;

% % Normalize original data
% MRNA_ORIGINAL = normalize(MRNA_ORIGINAL);
% PROTEIN_ORIGINAL = normalize(PROTEIN_ORIGINAL);
end




