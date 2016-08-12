function drawCheckDataDistribution(data,idx,text)
colors = {'g','b','r','c','m','y'};
marks = {'o','*','+'};
[~,~,V] = svd(data);
projection = data*V;
figure, clf
hold on
for i = 1:length(data)
%     plot(data(i,1),data(i,2),strcat(colors{mod(idx(i),6)+1},marks{ceil(idx(i)/6)}));
    plot(projection(i,1),projection(i,2),strcat(colors{mod(idx(i),6)+1},marks{ceil(idx(i)/6)}));
%     plot3(projection(i,1),projection(i,2),projection(i,3),strcat(colors{mod(idx(i),6)+1},marks{ceil(idx(i)/6)}));
end
title(text, 'FontSize', 30)
% xlabel('component1', 'FontSize', 20);
% ylabel('component2', 'FontSize', 20);
hold off
end