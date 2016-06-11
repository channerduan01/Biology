%
%
%

clear
load hmec;

% T X 6 arrays, M and P
%
mrna = data{1};
P1 = data{2}; protein = P1(:,2:7);
N = size(mrna,1);

figure(1), clf
subplot(121), imagesc([mrna zeros(N,1) protein]), colorbar('vert');
drawnow; hold on

c = zeros(N,1);
for n = 1:N
  cc = corrcoef(mrna(n,:), protein(n,:));
  c(n) = cc(1,2);
end

[xx,ii] = sort(c);
subplot(122), imagesc([mrna(ii,:) zeros(N,1) protein(ii,:)]);
colorbar('vert');
drawnow
%%
figure(2), clf,
bar(xx);
title('Correlations: mRNA and Protein', 'FontSize', 16);
grid on
%%
figure(3), clf,
plot(1:6,mrna(ii(1),:),'r', 1:6,protein(ii(1),:),'r-.',...
     1:6,mrna(ii(N),:),'b', 1:6,protein(ii(N),:),'b-.',...
     1:6,mrna(ii(N/2),:),'g', 1:6,protein(ii(N/2),:),'g-.',...
     'LineWidth',2);
grid on
title('Least (red) and Most (blue) Correlated Profiles', 'Fontsize', 16);
%%
l75 = length(find(c>0.75));

figure(1), clf,
imagesc([mrna(ii(N-l75+1:N),:) zeros(l75,1) protein(ii(N-l75+1:N),:)]);
colorbar('vert');
drawnow

for n=1:l75
    disp(names.GeneName{ii(N-n+1)})
end

    
    
