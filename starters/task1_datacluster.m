% Task 1, compare k-means and matrix approximate.
%
% init data
close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));


scene = 1;
switch scene
    case 1 % this case cannot be separate by NMF, but easy for k-means
        m1=[4;4];
        m2=[11,11];
        C1=[1,0;0,1];
        C2=[1,0;0,1];
    case 2 % this case can be separete by both
        m1=[4;12];
        m2=[12,4];
        C1=[1,0;0,1];
        C2=[1,0;0,1];
    case 3 % over the rank, NMF fails
        m1=[4;12];
        m2=[12,4];
        m3=[8;8];
        C1=[1,0;0,1];
        C2=[1,0;0,1];
        C3=[1,0;0,1];
    case 4 % more complex case, separable data
        m1=[4;4];
        m2=[8;7];
        m3=[5;10];
        C1=[1,0;0,1];
        C2=[2,1;1,2];
        C3=[2,-1;-1,2];
    case 5 % more complex case, inseparable data
        m1=[5;5];
        m2=[6;6];
        m3=[6;6];    
        C1=[1,0;0,1];
        C2=[2,1;1,2];
        C3=[2,-1;-1,2];
    case 6 % extremely seperable data! NMF wins!
        m1=[0,1];
        m2=[2,0];
        C1=[0,0;0,1];
        C2=[3,0;0,0];
    case 7 % NMF fail
        m1=[0,1];
        m2=[0,0];
        C1=[0,0;0,1];
        C2=[3,0;0,0];            
end
k = 2;
N = 80;
X1 = abs(mvnrnd(m1, C1, N));
X2 = abs(mvnrnd(m2, C2, N));
X1 = X1 + randn(size(X1))*0.1;    % adding noise
X2 = X2 + randn(size(X2))*0.1;    % adding noise

X = [X1;X2];
% NMF is sensitive to normalization!
% X = normalize(X);
% X = normalize_v2(X')';

X1 = X(1:N,:);
X2 = X(N+1:2*N,:);

% X1 = abs(X1);
% X2 = abs(X2);
if exist('m3','var'), X3 = mvnrnd(m3, C3, N); end;
figure(1), clf,
% grid on
hold on
plot(X1(:,1),X1(:,2),'bx');
plot(X2(:,1),X2(:,2),'ro');
if exist('X3','var'), plot(X3(:,1),X3(:,2),'g+'); end;
title('Original data', 'FontSize', 20)
xlabel('feature1', 'FontSize', 20);
ylabel('feature2', 'FontSize', 20);
set(gca,'FontSize',20);
if exist('X3','var')
    legend('class1', 'class2', 'class3');
else
    legend('class1', 'class2');
end
hold off

if exist('X3','var')
    data = [X1;X2;X3];
else
    data = [X1;X2];
end
ii = randperm(length(data));
X = data(ii,:);
% X = normalize(data(ii,:));
A = X';
% figure();
% plot(X(:,1),X(:,2),'mx');

cost = @(A,W,H) norm(A-W*H,'fro');
sparsity = @(L1,L2,dim) (sqrt(dim)-L1/sqrt(L2))/(sqrt(dim)-1);


%% k-means
repeat = 1;
opts = statset('Display','final');
[idx,centres] = kmeans(A',k,'Replicates',repeat,'Options',opts,'Start','sample');
W = centres';
H = zeros(k,length(idx));
for i = 1:length(idx)
    H(idx(i),i) = 1;
end
drawClusters(idx,centres,X,N,ii,'K-means');
fprintf('k-means-cost: %f\n', cost(A,W,H));
%% MU
% % [W,H,~] = mynmf(A,k,'verbose',0,'ALPHA',0,'BETA',0,'W_INIT',eye(k),'H_INIT',A);
% [W,H,~] = mynmf(A,k,'verbose',0,'ALPHA',2,'BETA',2);
% [~,idx] = max(H);
% drawClusters(idx',W',X,N,ii,'MU');
% fprintf('MU-cost: %f\n', cost(A,W,H));
% %% ALS
% [W,H,~] = mynmf(A,k,'METHOD','ALS','verbose',0,'ALPHA',1,'BETA',1);
% [~,idx] = max(H);
% drawClusters(idx',W',X,N,ii,'ALS');
% fprintf('ALS-cost: %f\n', cost(A,W,H));
% %% SNMF
% [W,H,iter,HIS] = nmf(A,k,'type','regularized','nnls_solver','bp','verbose',0,'ALPHA',1,'BETA',1);
% [~,idx] = max(H);
% % disp([norm(W,'fro'),norm(H,'fro')]);
% drawClusters(idx',W',X,N,ii,'NMF(BP)');
% fprintf('BP-cost: %f\n', cost(A,W,H));
%% NMFSC
[W,H,~] = mynmf(A,k,'METHOD','NMFSC','verbose',1,'ALPHA',1,'BETA',1,'RATE',0.5, 'MAX_ITER',80);
[~,idx] = max(H);
drawClusters(idx',W',X,N,ii,'NMFSC')
cost(A,W,H)

%% SVD ,,, cool!
% % information is reducing quicker and quicker while giving energy up
% A = randn(1000,1000);
% [U,S,V] = svd(A);
% norm(A-U*S*V','fro')
% for i = 1000:-1:50
%     S(i,:) = 0;
%     if mod(i,20) == 0
%         norm(U*S*V'-A,'fro')
%         rank(U*S*V')
%     end
% end
% for i = 50:-1:1
%     S(i,:) = 0;
%     norm(U*S*V'-A,'fro')
%     rank(U*S*V')
% end

%%
% [U,S,V] = svd(A);
% % S(:,3:400)=0;
% % V(:,3:400)=0;
% norm(A-U*S*V','fro')

%% Consistensy Analysis
% consistensyAnalysis(A,2:2,100,@wrapKmeanAsNmf)
% consistensyAnalysis(A,2:2,100,@mynmf)
% consistensyAnalysis(A,2:2,100,@(A,k) mynmf(A,k,'METHOD','ALS','MAX_ITER',100))
% consistensyAnalysis(A,2:2,100,@(A,k) nmf(A,k,'type','sparse','nnls_solver','bp','BETA',0))


%% Sparsity test
% L1 = 1;
% L2 = 0.6;
% sparsity(L1,L2,2)

% projection_operator(H(:,1),L1,L2)
% profile viewer







