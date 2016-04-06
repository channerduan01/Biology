


%% Draw figures of consistency

load consis

figure,clf
P = P_2_100_MRNA;
plot(2:100,P');
title('Consistency of mRNA', 'FontSize', 20)
xlabel('k', 'FontSize', 20);
ylabel('p', 'FontSize', 20);
set(gca,'FontSize',16);
legend('K-means', 'SNMF', 'MU', 'ALS', 'ALS-W');



figure,clf
P = P_2_100_CONC;
plot(2:100,P');
title('Consistency of concatenated mRNA and Protein', 'FontSize', 20)
xlabel('k', 'FontSize', 20);
ylabel('p', 'FontSize', 20);
set(gca,'FontSize',16);
legend('K-means', 'SNMF', 'MU', 'ALS', 'ALS-W');


