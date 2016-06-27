close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

K = 5;
J = K;

T = 6;
N = 500;

%%
AVG_K_ = rand(T, K);
COV_K_ = zeros(T, T, K);
for k = 1:K
%     COV_K_(:,:,k) = eye(T);
       COV_K_(:,:,k) = eye(T)*rand();     % easy distribution
%        COV_K_(:,:,k) = diag(rand(T,1));   % harder distribution
%     COV_K_(:,:,k) = rand(T,T);   % much harder distribution
%     COV_K_(:,:,k) = COV_K_(:,:,k) * COV_K_(:,:,k)';
    
end

AVG_J_ = rand(T,J);
COV_J_ = zeros(T, T, J);
for j = 1:J
%     COV_J_(:,:,j) = eye(T);
       COV_J_(:,:,j) = eye(T)*rand();
%        COV_J_(:,:,j) = diag(rand(T,1));
%     COV_J_(:,:,j) = rand(T,T);
%     COV_J_(:,:,j) = COV_J_(:,:,j) * COV_J_(:,:,j)';
end

% force pattern ==================
factor_gap = 1;
for k = 1:K
    AVG_K_(:,k) = ones(T,1)*k*factor_gap;
    AVG_J_(:,k) = ones(T,1)*k*factor_gap;
end
% ================================



MRNA = zeros(T, N);
PROTEIN = zeros(T, N);
idx_last = 1;
for k = 1:K
    num = floor((N-idx_last+1)/(K-k+1));
    range = idx_last:idx_last+num-1;
    MRNA(:, range) = mvnrnd(AVG_K_(:,k)', COV_K_(:,:,k), num)';
    PROTEIN(:, range) = mvnrnd(AVG_J_(:,k)', COV_J_(:,:,k), num)';
    idx_last = idx_last+num;
end

% PROTEIN = MRNA;

% permuted !!! important!!! ===================
% ii = randperm(N);
% PROTEIN = PROTEIN(:,ii);
% =============================================


% Calculate relationship between mRNA and protein
H1_ORIGINAL = zeros(K, N);
H2_ORIGINAL = zeros(J, N);
for k = 1:K
    H1_ORIGINAL(k, :) = mvnpdf(MRNA',AVG_K_(:,k)', COV_K_(:,:,k));
    H2_ORIGINAL(k, :) = mvnpdf(PROTEIN',AVG_J_(:,k)', COV_J_(:,:,k));
end
for i = 1:N
    H1_ORIGINAL(:, i) = H1_ORIGINAL(:, i)/sum(H1_ORIGINAL(:, i));
    H2_ORIGINAL(:, i) = H2_ORIGINAL(:, i)/sum(H2_ORIGINAL(:, i));
end

THETA_ORIGINAL = CalcuTheta(H1_ORIGINAL, H2_ORIGINAL, K, J, N);

figure();
imagesc(THETA_ORIGINAL);
title('Real THETA');


% Add noise
MRNA = MRNA + randn(size(MRNA))*1.8;
PROTEIN = PROTEIN + randn(size(PROTEIN))*1.5;

% Normalize original data
MRNA = normalize(MRNA);
PROTEIN = normalize(PROTEIN);
PROTEIN_ORIGINAL = PROTEIN;









