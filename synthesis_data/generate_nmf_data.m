%%
% Using the method from paper:
% Sparse Nonnegative Matrix Factorization for Clustering
%

close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

K = 6;
J = K;

T = 100;
N = 1000;

VALUE_RANGE_D = 1:3;
VALUE_VARIANCE = 0.5;


%%
AVG_K_ = zeros(T, K);
COV_K_ = zeros(T, T, K);
for extra = 1:1
    for i = 1:T
        k = floor(rand()*K)+1;
        AVG_K_(i,k) = VALUE_RANGE_D(floor(rand()*length(VALUE_RANGE_D))+1);
        COV_K_(i,i,k) = VALUE_VARIANCE;
    end
end
% for k = 1:K
%     for i = 1:100
%         t = floor(rand()*T)+1;
%         AVG_K_(t,k) = VALUE_RANGE_D(floor(rand()*length(VALUE_RANGE_D))+1);
%         COV_K_(t,t,k) = 0.3;
%     end
% end

% AVG_J_ = AVG_K_;
% COV_J_ = COV_K_;
AVG_J_ = zeros(T, J);
COV_J_ = zeros(T, T, J);
for extra = 1:1
    for i = 1:T
        j = floor(rand()*J)+1;
        AVG_J_(i,j) = VALUE_RANGE_D(floor(rand()*length(VALUE_RANGE_D))+1);
        COV_J_(i,i,j) = VALUE_VARIANCE;
    end
end

MRNA = zeros(T, N);
PROTEIN = zeros(T, N);
J_RANGE = randperm(K);
range = @(K,N,k) floor(N/K*(k-1))+1:floor(N/K*k);
for k = 1:K
    k_ = J_RANGE(k);
    MRNA(:, range(K,N,k)) = mvnrnd(AVG_K_(:,k)', COV_K_(:,:,k), length(range(K,N,k)))';
    PROTEIN(:, range(K,N,k)) = mvnrnd(AVG_J_(:,k_)', COV_J_(:,:,k_), length(range(K,N,k)))';    
end

% Calculate relationship between mRNA and protein
H1_ORIGINAL = zeros(K, N);
H2_ORIGINAL = zeros(J, N);
for k = 1:K
    H1_ORIGINAL(k, :) = MyMvnpdf(MRNA',AVG_K_(:,k)', COV_K_(:,:,k));
    H2_ORIGINAL(k, :) = MyMvnpdf(PROTEIN',AVG_J_(:,k)', COV_J_(:,:,k));
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
% PROTEIN = PROTEIN(:,ii);
% =============================================


% Add noise =============================
% MRNA = MRNA + randn(size(MRNA))*2;
% PROTEIN = PROTEIN + randn(size(PROTEIN))*2;

% MRNA = MRNA + randn(size(MRNA))*0.1;
% PROTEIN = PROTEIN + randn(size(PROTEIN))*1;

% MRNA = MRNA + randn(size(MRNA))*2;
PROTEIN = PROTEIN + randn(size(PROTEIN))*2;

MRNA = MRNA + rand(size(MRNA))*2;
PROTEIN = PROTEIN + rand(size(PROTEIN))*6;

% MRNA = MRNA + wgn(size(MRNA,1),size(MRNA,2),1)*1;
% PROTEIN = PROTEIN + wgn(size(PROTEIN,1),size(PROTEIN,2),1)*1;

THETA_ORIGINAL = CalcuTheta(H1_ORIGINAL, H2_ORIGINAL, K, J, N);
% =======================================


% MRNA = MRNA-3;

% Normalize original data
% MRNA = normalize(MRNA);
% PROTEIN = normalize(PROTEIN);
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



