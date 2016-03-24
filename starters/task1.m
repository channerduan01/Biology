% Task 1, compare k-means and matrix approximate.
%

%% init data
close all
clear
clc
% separable data
% m1=[4;4];
% m2=[8;7];
% m3=[4;10];
% inseparable data
m1=[5;5];
m2=[6;6];
m3=[6;6];

C1=[1,0;0,1];
C2=[2,1;1,2];
C3=[2,-1;-1,2];

N = 500;
X1 = mvnrnd(m1, C1, N);
X2 = mvnrnd(m2, C2, N);
X3 = mvnrnd(m3, C3, N);
figure(1), clf,
% grid on
hold on
plot(X1(:,1),X1(:,2),'bx');
plot(X2(:,1),X2(:,2),'ro');
plot(X3(:,1),X3(:,2),'g+');
title('Original data', 'FontSize', 20)
xlabel('feature1', 'FontSize', 20);
ylabel('feature2', 'FontSize', 20);
axis([0 15 0 15]);
set(gca,'FontSize',20);
legend('class1', 'class2', 'class3');
hold off
data = [X1;X2;X3];

%% k-means
ii = randperm(length(data));
X = data(ii,:);
k = 3;
opts = statset('Display','final');
[idx,centres] = kmeans(X,k,'Replicates',20,'Options',opts,'Start','sample');
figure, clf,
hold on
for i = 1:length(idx)
    str = '';
    switch idx(i)
        case 1
            str = strcat(str,'b');
        case 2
            str = strcat(str,'r');
        case 3            
            str = strcat(str,'g');
    end
    switch floor((ii(i)-1)/N)
        case 0
            str = strcat(str,'x');
        case 1
            str = strcat(str,'o');
        case 2
            str = strcat(str,'+');
    end
    plot(X(i,1),X(i,2),str);
end
plot(centres(:,1),centres(:,2),'kx','MarkerSize',15,'LineWidth',5);
title('k-means', 'FontSize', 20)
xlabel('feature1', 'FontSize', 20);
ylabel('feature2', 'FontSize', 20);
axis([0 15 0 15]);
hold off






