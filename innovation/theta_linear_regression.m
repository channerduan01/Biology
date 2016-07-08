%% research for calculate Theta by linear regression

close all;
clear;
clc;

% init
K = 15;
J = 19;
T = 300;

W1 = rand(T,K);
W2 = rand(T,J);
%%
THETA1 = rand(K,J);
cost = @(V,W,H) sqrt(sum(sum((V-W*H).^2)));


fprintf('init cost: %f\n', cost(W2,W1,THETA1));
alpha = 0.0001;
for i = 1:50
    THETA1 = 0.9*THETA1 - alpha * W1'*(W1*THETA1-W2);
    fprintf('iter-%d cost: %f\n', i, cost(W2,W1,THETA1));
end




