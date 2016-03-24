% Task 1, compare k-means and matrix approximate.
%

%% init data
close all
clear
clc

m1=[4;4];
m2=[8;7];
m3=[4;10];
C1=[1,0;0,1];
C2=[2,1;1,2];
C3=[2,-1;-1,2];

N = 200;
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

%% k-means

