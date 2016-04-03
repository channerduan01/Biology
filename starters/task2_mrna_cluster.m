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

A = mrna';
% A = conca';

min_ = min(min(A));
max_ = max(max(A));
A = (A-min_)/(max_-min_);
% figure(1), clf, imagesc([mrna zeros(N,1) A'*(max_-min_)+min_])

k = 9;
cost = @(A,W,H) sqrt(sum(sum((A-W*H).^2)));

%% k-means
repeat = 1;
opts = statset('Display','final');
[idx,centres] = kmeans(A',k,'Replicates',repeat,'Options',opts,'Start','sample');
W = centres';
H = zeros(k,length(idx));
for i = 1:length(idx)
    H(idx(i),i) = 1;
end
disp(['k-means cost:',num2str(cost(A,W,H))]);
drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'K-means');
%% MU
[W,H] = mynmf(A,k,'verbose',1,'MAX_ITER',100);
[~,idx] = max(H);
cost(A,W,H)
drawGeneTimesequence(A',idx,k);
% drawCheckDataDistribution(A',idx,'MU');
%% ALS
[W,H] = mynmf(A,k,'METHOD','ALS','verbose',1);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'ALS');
%% ANLS
[W,H] = nmf(A,k,'type','sparse','nnls_solver','bp','verbose',1);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'ANLS');

%% Consistency Analysis
repeatTime = 2;
K = 10:13;
P = zeros(3,length(K));
P(1,:) = consistensyAnalysis(A,K,repeatTime,@wrapKmeanAsNmf);
P(2,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) nmf(A,k,'type','sparse','MAX_ITER',100));
P(3,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) mynmf(A,k,'METHOD','MU','MAX_ITER',100));
P(4,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) mynmf(A,k,'METHOD','ALS','MAX_ITER',100));
P(5,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) mynmf(A,k,'METHOD','ALS_W','MAX_ITER',100));

% Draw figures of consistency
figure,clf
plot(K,P');
title('Consistency', 'FontSize', 20)
xlabel('k', 'FontSize', 20);
ylabel('p', 'FontSize', 20);
set(gca,'FontSize',16);
legend('K-means', 'Sparse', 'MU', 'ALS', 'ALS-W');




