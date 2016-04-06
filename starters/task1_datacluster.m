% Task 1, compare k-means and matrix approximate.
%
% init data
close all
clear
clc

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
end
k = 2;
N = 200;
X1 = mvnrnd(m1, C1, N);
X2 = mvnrnd(m2, C2, N);
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
axis([0 15 0 15]);
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
A = X';

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
drawClusters(idx,centres,X,N,ii,'K-means')

%% MU
[W,H,~] = mynmf(A,k,'verbose',0);
[~,idx] = max(H);
drawClusters(idx',W',X,N,ii,'MU')
cost(A,W,H)
%% ALS
[W,H,~] = mynmf(A,k,'METHOD','ALS','verbose',0);
[~,idx] = max(H);
drawClusters(idx',W',X,N,ii,'ALS')
cost(A,W,H)
%% CVX
% [W,H,~] = mynmf(A,k,'METHOD','CVX','verbose',1,'MAX_ITER',20);
% [~,idx] = max(H);
% drawClusters(idx',W',X,N,ii,'CVX')
% cost(A,W,H)
%% SNMF
[W,H,iter,HIS] = nmf(A,k,'type','sparse','nnls_solver','bp');
[~,idx] = max(H);
drawClusters(idx',W',X,N,ii,'SNMF')
cost(A,W,H)

%% SVD ,,, cool!
A = randn(1000,1000);
[U,S,V] = svd(A);
norm(A-U*S*V','fro')
for i = 1000:-1:50
    S(i,:) = 0;
    if mod(i,20) == 0
        norm(U*S*V'-A,'fro')
        rank(U*S*V')
    end
end
for i = 50:-1:1
    S(i,:) = 0;
    norm(U*S*V'-A,'fro')
    rank(U*S*V')
end

%%
[U,S,V] = svd(A);
S(:,3:400)=0;
% V(:,3:400)=0;
norm(A-U*S*V','fro')

%% Consistensy Analysis
% consistensyAnalysis(A,2:2,100,@wrapKmeanAsNmf)
consistensyAnalysis(A,2:2,100,@mynmf)
consistensyAnalysis(A,2:2,100,@(A,k) mynmf(A,k,'METHOD','ALS','MAX_ITER',100))
% consistensyAnalysis(A,2:2,100,@(A,k) nmf(A,k,'type','sparse','nnls_solver','bp','BETA',0))










