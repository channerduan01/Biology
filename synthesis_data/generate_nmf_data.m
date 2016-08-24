%%
% Using the method from paper:
% Sparse Nonnegative Matrix Factorization for Clustering
%

close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

K = 5;
J = K;

T = 10;
N = 1000;

VALUE_RANGE_D = 1:0.1:2;
VALUE_VARIANCE = 0.3;


%%
while true
    AVG_K_ = zeros(T, K);
    COV_K_ = zeros(T, T, K);
    for extra = 1:1
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
    for extra = 1:1
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

    

MRNA_ORIGINAL = zeros(T, N);
PROTEIN_ORIGINAL = zeros(T, N);
J_RANGE = randperm(K);
range = @(K,N,k) floor(N/K*(k-1))+1:floor(N/K*k);
for k = 1:K
    k_ = J_RANGE(k);
    MRNA_ORIGINAL(:, range(K,N,k)) = mvnrnd(AVG_K_(:,k)', COV_K_(:,:,k), length(range(K,N,k)))';
    PROTEIN_ORIGINAL(:, range(K,N,k)) = mvnrnd(AVG_J_(:,k_)', COV_J_(:,:,k_), length(range(K,N,k)))';    
end


MRNA_ORIGINAL(MRNA_ORIGINAL<0) = 0;
PROTEIN_ORIGINAL(PROTEIN_ORIGINAL<0) = 0;

% Calculate relationship between MRNA_ORIGINAL and PROTEIN_ORIGINAL
H1_ORIGINAL = zeros(K, N);
H2_ORIGINAL = zeros(J, N);
for k = 1:K
    H1_ORIGINAL(k, :) = MyMvnpdf(MRNA_ORIGINAL',AVG_K_(:,k)', COV_K_(:,:,k));
    H2_ORIGINAL(k, :) = MyMvnpdf(PROTEIN_ORIGINAL',AVG_J_(:,k)', COV_J_(:,:,k));
end
for i = 1:N
    H1_ORIGINAL(:, i) = H1_ORIGINAL(:, i)/sum(H1_ORIGINAL(:, i));
    H2_ORIGINAL(:, i) = H2_ORIGINAL(:, i)/sum(H2_ORIGINAL(:, i));
end

THETA_ORIGINAL = CalcuTheta(H1_ORIGINAL, H2_ORIGINAL, K, J, N);

figure();
imagesc(THETA_ORIGINAL);
title('Real THETA');


% permuted !!! important!!! ===================
% ii = randperm(N);
% PROTEIN_ORIGINAL = PROTEIN_ORIGINAL(:,ii);
% =============================================

%%
MRNA = MRNA_ORIGINAL;
PROTEIN = PROTEIN_ORIGINAL;

% Add noise                ==============
amplitude = 1;
NOISE_1 = randn(size(MRNA))*amplitude;
NOISE_2 = randn(size(PROTEIN))*amplitude;
MRNA(MRNA>0) = MRNA(MRNA>0) + NOISE_1(MRNA>0);
PROTEIN(PROTEIN>0) = PROTEIN(PROTEIN>0) + NOISE_2(PROTEIN>0);

% MRNA = MRNA + rand(size(MRNA))*1;
% PROTEIN = PROTEIN + rand(size(PROTEIN))*1;
% =======================================


% Normalize original data
MRNA = normalize(MRNA);
PROTEIN = normalize(PROTEIN);

% MRNA = normalize_v2(MRNA);
% PROTEIN = normalize_v2(PROTEIN);
% PROTEIN_ORIGINAL = PROTEIN;
% MRNA(MRNA<0) = 0;
% PROTEIN(PROTEIN<0) = 0;
% 
% MRNA = pow2(MRNA);
% PROTEIN = pow2(PROTEIN);
% MRNA(MRNA>0.5) = 0;
% PROTEIN(PROTEIN>0.5) = 0;

%% visualize the distribution
[~,idx_] = max(H1_ORIGINAL);
drawCheckDataDistribution(MRNA',idx_,'MRNA distribution');





%%
C_ = CalcuCovariance(MRNA');
sum(sum(abs(C_)))









