function drawGeneTimesequence(data,idx,k)
figure, clf
set(gca,'FontSize',20);
mode = 100 + k*10;
for i = 1:k
    subplot(mode+i), imagesc(data(idx==i,:)), title(['Num:',num2str(sum(idx==i))])
end
set(gcf,'position',[100 50 1300 750])
drawnow;
end