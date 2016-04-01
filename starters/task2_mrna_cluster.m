% Task 1, compare k-means and matrix approximate.
%
% init data
close all
clear
clc

load hmec;
mrna = data{1};
P1 = data{2}; 
protein = P1(:,2:7);
conca = [mrna,protein];
N = size(mrna,1);

% A = mrna';
A = conca';
min_ = min(min(A));
max_ = max(max(A));
A = (A-min_)/(max_-min_);
% figure(1), clf, imagesc([mrna zeros(N,1) A'*(max_-min_)+min_])

k = 15;
cost = @(A,W,H) sqrt(sum(sum((A-W*H).^2)));

%% k-means
repeat = 10;
opts = statset('Display','final');
[idx,centres] = kmeans(A',k,'Replicates',repeat,'Options',opts,'Start','sample');
W = centres';
H = zeros(k,length(idx));
for i = 1:length(idx)
    H(idx(i),i) = 1;
end
disp(['k-means cost:',num2str(cost(A,W,H))]);
% drawGeneTimesequence(mrna,idx,k);
drawCheckDataDistribution(mrna,idx,'K-means');
%% MU
[W,H,~,HIS] = mynmf(A,k,'verbose',1,'MAX_ITER',100);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(mrna,idx,k);
drawCheckDataDistribution(mrna,idx,'MU');
%% ALS
[W,H,~] = mynmf(A,k,'METHOD','ALS','verbose',1);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(mrna,idx,k);
drawCheckDataDistribution(mrna,idx,'ALS');
%% ANLS
[W,H,iter,HIS] = nmf(A,k,'type','sparse','nnls_solver','bp','verbose',1);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(mrna,idx,k);
drawCheckDataDistribution(mrna,idx,'ANLS');















