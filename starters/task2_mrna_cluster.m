% Task 2
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

% min_ = min(min(A));
% max_ = max(max(A));
% A = (A-min_)/(max_-min_);
A = normalize(A);
% figure(1), clf, imagesc([mrna zeros(N,1) A'*(max_-min_)+min_])

k = 2;
cost = @(A,W,H) norm(A-W*H,'fro');

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
% drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'K-means');
% cost(A,W,H)
%% MU
[W,H] = mynmf(A,k,'verbose',0,'MAX_ITER',100);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'MU');
%% ALS
[W,H] = mynmf(A,k,'METHOD','ALS','ALPHA',0.1,'BETA',0.1,'verbose',0,'ALPHA',1,'BETA',1);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'ALS');
% %% ALS-W
% [W,H] = mynmf(A,k,'METHOD','ALS_W','verbose',0);
% [~,idx] = max(H);
% cost(A,W,H)
% % drawGeneTimesequence(A',idx,k);
% drawCheckDataDistribution(A',idx,'ALS-W');
%% SNMF
[W,H] = nmf(A,k,'type','sparse','ALPHA',0.1,'BETA',0.1,'nnls_solver','bp','MAX_ITER',100,'verbose',1);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'SNMF');
%% NMFSC
[W,H,~] = mynmf(A,k,'METHOD','NMFSC','verbose',1,'ALPHA',1,'BETA',1,'RATE',10,'MAX_ITER',5);
[~,idx] = max(H);
cost(A,W,H)
% drawGeneTimesequence(A',idx,k);
drawCheckDataDistribution(A',idx,'NMFSC');

%% SVD ,,, confusing part~
[U,S,V] = svd(A);
S(k+1:size(S,1),:) = 0;
norm(A-U*S*V','fro')


%% Consistency Analysis
tic;
repeatTime = 100;
K = 2:2;
P = zeros(6,length(K));
P(1,:) = consistensyAnalysis(A,K,repeatTime,@wrapKmeanAsNmf);
P(2,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) nmf(A,k,'type','sparse','MAX_ITER',100));
% P(3,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) mynmf(A,k,'METHOD','MU','MAX_ITER',100));
% P(4,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) mynmf(A,k,'METHOD','ALS','MAX_ITER',100,'ALPHA',0,'BETA',10));
% P(5,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) mynmf(A,k,'METHOD','ALS_W','MAX_ITER',100));
% P(6,:) = consistensyAnalysis(A,K,repeatTime,@(A,k) mynmf(A,k,'METHOD','NMFSC','ALPHA',1,'BETA',1,'RATE',0.5,'MAX_ITER',20));
toc;

%%
% Draw figures of consistency
% figure,clf
% P = P_2_100_MRNA;
% plot(2:100,P');
% title('Consistency', 'FontSize', 20)
% xlabel('k', 'FontSize', 20);
% ylabel('p', 'FontSize', 20);
% set(gca,'FontSize',16);
% legend('K-means', 'SNMF', 'MU', 'ALS', 'ALS-W');


%%

