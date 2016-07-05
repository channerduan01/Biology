function drawClusters(idx,centres,X,N,ii,text)
colors = {'b','r','g','c','m','y'};
original_marks = {'x','o','+'};
figure, clf,
hold on
for i = 1:length(idx)
    plot(X(i,1),X(i,2),strcat(colors{idx(i)},original_marks{floor((ii(i)-1)/N)+1}));
end
plot(centres(:,1),centres(:,2),'kx','MarkerSize',15,'LineWidth',5);
title(text, 'FontSize', 20)
xlabel('feature1', 'FontSize', 20);
ylabel('feature2', 'FontSize', 20);
% legend('class1', 'class2');
% axis([0 15 0 15]);
hold off
end