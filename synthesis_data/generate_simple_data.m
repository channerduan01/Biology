close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

K = 5;
J = K;

T = 10;
N = 1000;

%%
% AVG_K_ = rand(T, K);
% AVG_J_ = rand(T, J);
gap_requirement = 0.07;

while true
    AVG_K_ = rand(T, K);
    pass_ = true;
    for t = 2:T
        for tt = 1:t-1
            if norm(AVG_K_(tt, K)-AVG_K_(t, K), 'fro') < gap_requirement
                pass_ = false;
                break;
            end
        end
        if ~pass_, break; end
    end
    if pass_, break; end
end
while true
    AVG_J_ = rand(T, J);
    t = 1;
    pass_ = true;
    for t = 2:T
        for tt = 1:t-1
            if  norm(AVG_J_(tt, J)-AVG_J_(t, J), 'fro') < gap_requirement
                pass_ = false;
                break;
            end
        end
        if ~pass_, break; end
    end
    if pass_, break; end
end

% t = 1;
% while true
%     AVG_K_(t, :) = rand(1, K);
%     AVG_J_(t, :) = rand(1, J);
%     pass_ = true;
%     for tt = 1:t-1
%         if norm(AVG_K_(tt, K)-AVG_K_(t, K), 'fro') < gap_requirement || ...
%             norm(AVG_J_(tt, J)-AVG_J_(t, J), 'fro') < gap_requirement
%             pass_ = false;
%             break;
%         end
%     end
%     if pass_
%         t = t+1;
%         if t > T, break; end
%     end
% end

% force pattern ==================
% factor_gap = 1;
% for k = 1:K
%     AVG_K_(:,k) = ones(T,1)*k*factor_gap;
%     AVG_J_(:,k) = ones(T,1)*k*factor_gap;
% end
% ================================

COV_K_ = zeros(T, T, K);
for k = 1:K
%     COV_K_(:,:,k) = eye(T);
%        COV_K_(:,:,k) = eye(T)*rand();     % easy distribution
       COV_K_(:,:,k) = diag(rand(T,1))/2;   % harder distribution
%     COV_K_(:,:,k) = rand(T,T);   % much harder distribution
%     COV_K_(:,:,k) = COV_K_(:,:,k) * COV_K_(:,:,k)';
    
end

COV_J_ = zeros(T, T, J);
for j = 1:J
%     COV_J_(:,:,j) = eye(T);
%        COV_J_(:,:,j) = eye(T)*rand();
       COV_J_(:,:,j) = diag(rand(T,1))/2;
%     COV_J_(:,:,j) = rand(T,T);
%     COV_J_(:,:,j) = COV_J_(:,:,j) * COV_J_(:,:,j)';
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

% for nonnegative
% MRNA(MRNA<0) = 0;
% PROTEIN(PROTEIN<0) = 0;

% idx_last = 1;
% for k = 1:K
%     num = floor((N-idx_last+1)/(K-k+1));
%     range = idx_last:idx_last+num-1;
%     MRNA(:, range) = mvnrnd(AVG_K_(:,k)', COV_K_(:,:,k), num)';
%     PROTEIN(:, range) = mvnrnd(AVG_J_(:,k)', COV_J_(:,:,k), num)';
%     idx_last = idx_last+num;
% end

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


% Add noise =============================
% MRNA = MRNA + randn(size(MRNA))*1;
% PROTEIN = PROTEIN + randn(size(PROTEIN))*0.1;

% MRNA = MRNA + randn(size(MRNA))*0.1;
% PROTEIN = PROTEIN + randn(size(PROTEIN))*1;

% MRNA = MRNA + randn(size(MRNA))*1;
% PROTEIN = PROTEIN + randn(size(PROTEIN))*1;
% =======================================


% Normalize original data
% MRNA = normalize(MRNA);
% PROTEIN = normalize(PROTEIN);


% PROTEIN_ORIGINAL = PROTEIN;
% MRNA = abs(MRNA);
% PROTEIN = abs(PROTEIN);





