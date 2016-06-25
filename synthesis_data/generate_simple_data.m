close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

K = 15;
J = 15;

T = 100;
N = 500;

%%
AVG_K = rand(T, K);
COV_K = zeros(T, T, K);
for k = 1:K
   COV_K(:,:,k) = eye(T);
end

AVG_J = rand(T,J);
COV_J = zeros(T, T, J);
for j = 1:J
   COV_K(:,:,j) = eye(T); 
end


MRNA = zeros(T, N);
PROTEIN = zeros(T, N);

H1_ORIGINAL = zeros(K, N);
H2_ORIGINAL = zeros(J, N);

idx_last = 1;
for k = 1:K
    num = floor((N-idx_last+1)/(K-k+1));
%     fprintf('k:%d idx_last:%d num:%d \n',k,idx_last,num);
    range = idx_last:idx_last+num-1;
    MRNA(:, range) = mvnrnd(AVG_K(:,k)', COV_K(:,:,k), num)';
    H1_ORIGINAL(k, range) = 1;
    H2_ORIGINAL(k, range) = 1;
    idx_last = idx_last+num;
end


PROTEIN = MRNA;
PROTEIN_ORIGINAL = MRNA;





