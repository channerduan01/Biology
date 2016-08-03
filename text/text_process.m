%% Roger EM coupled clustering model
%
% init data
close all
clc
clear
addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

load tfidf_train_data
V = tfidf(1:3000,:);
K = 5;

%%

[W,H,~] = mynmf(V',K,'verbose',1,'METHOD','ALS','ALPHA',1,'BETA',1,'MAX_ITER',3,'MIN_ITER',3);