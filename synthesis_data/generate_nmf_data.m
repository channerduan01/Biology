%%
% Using the method from paper:
% Sparse Nonnegative Matrix Factorization for Clustering
%

close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

K = 15;
J = K;

T = 300;
N = 1000;

VALUE_RANGE_D = 3:6;
VALUE_VARIANCE = 0.1;

%%
AVG_K_ = zeros(T, K);
COV_K_ = zeros(T, T, K);
for i = 1:T
    k = floor(rand()*K)+1;
%     k = i;
    AVG_K_(i,k) = VALUE_RANGE_D(floor(rand()*length(VALUE_RANGE_D))+1);
    COV_K_(i,i,k) = VALUE_VARIANCE;
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
for i = 1:T
    j = floor(rand()*J)+1;
    AVG_J_(i,j) = VALUE_RANGE_D(floor(rand()*length(VALUE_RANGE_D))+1);
    COV_J_(i,i,j) = VALUE_VARIANCE;
end

MRNA = zeros(T, N);
PROTEIN = zeros(T, N);
idx_last = 1;
for k = 1:K
    num = floor((N-idx_last+1)/(K-k+1));
    range = idx_last:idx_last+num-1;
%     MRNA(:, range) = abs(mvnrnd(AVG_K_(:,k)', COV_K_(:,:,k), num)');
%     PROTEIN(:, range) = abs(mvnrnd(AVG_J_(:,k)', COV_J_(:,:,k), num)');
    
    MRNA(:, range) = mvnrnd(AVG_K_(:,k)', COV_K_(:,:,k), num)';
    PROTEIN(:, range) = mvnrnd(AVG_J_(:,k)', COV_J_(:,:,k), num)';    
    
    idx_last = idx_last+num;
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




% Add noise =============================
MRNA = MRNA + randn(size(MRNA))*1;
PROTEIN = PROTEIN + randn(size(PROTEIN))*0.1;

% MRNA = MRNA + randn(size(MRNA))*0.1;
% PROTEIN = PROTEIN + randn(size(PROTEIN))*2;

% MRNA = MRNA + randn(size(MRNA))*1;
% PROTEIN = PROTEIN + randn(size(PROTEIN))*1;
% =======================================


% Normalize original data
% MRNA = normalize(MRNA);
% PROTEIN = normalize(PROTEIN);
% PROTEIN_ORIGINAL = PROTEIN;
% MRNA(MRNA<0) = 0;
% PROTEIN(PROTEIN<0) = 0;
% MRNA = pow2(MRNA);
% PROTEIN = pow2(PROTEIN);



