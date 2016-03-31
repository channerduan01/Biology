% Task 1, compare k-means and matrix approximate.
%

%% init data
close all
clear
clc


m1=[4;4];
m2=[11,11];
% m1=[4;12];
% m2=[12,4];
% m3=[8;8];
% separable data
% m1=[4;4];
% m2=[8;7];
% m3=[5;10];
% inseparable data
% m1=[5;5];
% m2=[6;6];
% m3=[6;6];

C1=[1,0;0,1];
% C2=[2,1;1,2];
C2=[1,0;0,1];
% C3=[2,-1;-1,2];
% C3=[1,0;0,1];

N = 200;
X1 = mvnrnd(m1, C1, N);
X2 = mvnrnd(m2, C2, N);
% X3 = mvnrnd(m3, C3, N);
figure(1), clf,
% grid on
hold on
plot(X1(:,1),X1(:,2),'bx');
plot(X2(:,1),X2(:,2),'ro');
% plot(X3(:,1),X3(:,2),'g+');
title('Original data', 'FontSize', 20)
xlabel('feature1', 'FontSize', 20);
ylabel('feature2', 'FontSize', 20);
axis([0 15 0 15]);
set(gca,'FontSize',20);
% legend('class1', 'class2', 'class3');
legend('class1', 'class2');
hold off

% data = [X1;X2;X3];
data = [X1;X2];
ii = randperm(length(data));
X = data(ii,:);
A = X';
k = 2;
r = k;

cost = @(A,W,H) sqrt(sum(sum((A-W*H).^2)));

%% k-means
repeat = 1;
opts = statset('Display','final');
[idx,centres] = kmeans(X,k,'Replicates',repeat,'Options',opts,'Start','sample');
W = centres';
H = zeros(k,length(idx));
for i = 1:length(idx)
    H(idx(i),i) = 1;
end

disp(['k-means cost:',num2str(cost(A,W,H))]);
drawClusters(idx,centres,X,N,ii,'k-means')

%% CNMF
a = mean(A(:));
b = mean(A(:));
[m,n] = size(A);

% Multiplicative update rules(Lee and Seung, 2001) and CNMF(Pauca et al., 2006)
W = rand(m,r);
H = rand(r,n);
max_iter = 50;
for i = 1:max_iter
    H = H .* (W'*A)./(W'*W*H + b*H + eps);
    %     H = H .* (W'*A)./(W'*W*H + b*ones(r)*H + eps);
    W = W .* (A*H')./(W*H*H' + a*W + eps);
end

[~,idx] = max(H);
idx = idx';
% centres = zeros(r,m);
% for i = 1:r
%     centres(i,:) = mean(X(idx==i,:));
% end
centres = W';
drawClusters(idx,centres,X,N,ii,'CNMF')

%% ALS
a = mean(A(:));
b = mean(A(:));
[m,n] = size(A);

% Basic ALS(Paatero and Tapper, 1994)
W = rand(m,r);
max_iter = 20;
for i = 1:max_iter
    H = pinv(W'*W+b*eye(r))*(W'*A);
    H(H<0)=0;
    W = (pinv(H*H'+a*eye(r))*(H*A'))';
    W(W<0)=0;
    %     if mod(i,10) == 0
    disp(['iter ',num2str(i),' - cost ',num2str(cost(A,W,H))]);
    %     end
end

[~,idx] = max(H);
idx = idx';
% centres = zeros(r,m);
% for i = 1:r
%     centres(i,:) = mean(X(idx==i,:));
% end
centres = W';
drawClusters(idx,centres,X,N,ii,'ALS')

%% ANLS
[W,H,iter,HIS] = nmf(A,k,'type','sparse','nnls_solver','bp','alpha',1.1,'beta',0.1);
[~,idx] = max(H);
drawClusters(idx',W',X,N,ii,'ANLS')






















