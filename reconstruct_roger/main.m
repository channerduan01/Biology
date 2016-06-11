% Roger EM coupled clustering model
%
% init data
close all
clear
clc

cost = @(V,W,H) norm(V-W*H,'fro');

load hmec;
mrna = data{1};
P1 = data{2}; 
protein = P1(:,2:7);
mrna = normalize(mrna)';
protein = normalize(protein');

%% K-mean test
% V = protein;
% K = 10;
% 
% [centres, ~] = mykmeans(V, K, 10);
% centres
% repeat = 1;
% opts = statset('Display','final');
% [idx,centres] = kmeans(V',K,'Replicates',repeat,'Options',opts,'Start','sample');
% W = centres';
% H = zeros(K,length(idx));
% for i = 1:length(idx)
%     H(idx(i),i) = 1;
% end
% disp(['k-means cost:',num2str(cost(V,W,H))]);

%% Roger's Model
K = 15;
J = 19;
MAX_ITER = 100;
patience = 2;
MRNA = mrna;
PROTEIN = protein;

for i = 1:20
    [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,low_bound] = ...
        mycoupleclustering(MRNA, PROTEIN, K, J, MAX_ITER, patience, false);
    low_bound
end





