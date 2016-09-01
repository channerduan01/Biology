%%
% Methods comparison

close all;
clear;
clc;
addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

T = 6;
N = 1000;
K = 4;
J = 4;


[MRNA_ORIGINAL, PROTEIN_ORIGINAL, THETA_ORIGINAL, H1_ORIGINAL, H2_ORIGINAL] = GenerateData(K, J, T, N);

% ttt = randperm(N);
% PROTEIN_ORIGINAL = PROTEIN_ORIGINAL(:, ttt);
% H2_ORIGINAL = H2_ORIGINAL(:, ttt);

%%
REPEAT_NUM = 20;
% set_of_amplitude = 0:0.2:2;
set_of_amplitude = [0,0];

HUGE_TALBE = zeros(REPEAT_NUM, 4*5, length(set_of_amplitude));
BEST_HUGE_TALBE = zeros(4, 3, length(set_of_amplitude));
ERROR_TABLE = zeros(4, 1, length(set_of_amplitude));
%% create specific noise
MRNA_ORIGINAL = MRNA;
PROTEIN_ORIGINAL = PROTEIN;
NOISE_1 = randn(size(MRNA_ORIGINAL));
NOISE_2 = randn(size(PROTEIN_ORIGINAL));


%%
for i = 1:length(set_of_amplitude)
fprintf('\n\n\n>>>>> deal with amplitude: %f\n', set_of_amplitude(i));

% % add noise
MRNA = MRNA_ORIGINAL;
PROTEIN = PROTEIN_ORIGINAL;
MRNA(MRNA>0) = MRNA(MRNA>0) + set_of_amplitude(i)*NOISE_1(MRNA>0);
PROTEIN(PROTEIN>0) = PROTEIN(PROTEIN>0) + set_of_amplitude(i)*NOISE_2(PROTEIN>0);
MRNA = normalize(MRNA);
PROTEIN = normalize(PROTEIN);

% % run algorithms
[HUGE_TALBE(:,:,i), BEST_HUGE_TALBE(:,:,i), ERROR_TABLE(:,:,i)] ...
    = RunAllAlgorithm(MRNA, PROTEIN, THETA_ORIGINAL, H1_ORIGINAL, H2_ORIGINAL, K, J, T, N, false ...
    , 5);
end


%% RESULT
result_ = reshape(mean(HUGE_TALBE), [20,length(set_of_amplitude)]);

figure;
hold on;
plot(set_of_amplitude, result_(2,:)*100);
plot(set_of_amplitude, result_(7,:)*100);
plot(set_of_amplitude, result_(12,:)*100);
plot(set_of_amplitude, result_(17,:)*100);
set(gca,'FontSize',20);
legend('Double K-means', 'Double NMF', 'Rogers'' Model', 'Coupled NMF');
hold off;
%% RESULT
result_ = reshape(mean(HUGE_TALBE), [20,length(set_of_amplitude)]);

figure;
hold on;
plot(set_of_amplitude, result_(3,:)*100);
plot(set_of_amplitude, result_(8,:)*100);
plot(set_of_amplitude, result_(13,:)*100);
plot(set_of_amplitude, result_(18,:)*100);
set(gca,'FontSize',20);
legend('Double K-means', 'Double NMF', 'Rogers'' Model', 'Coupled NMF');
hold off;
%% Best RESULT
best_result_ = reshape(BEST_HUGE_TALBE(:,1,:), [4,length(set_of_amplitude)]);

figure;
hold on;
plot(best_result_(1,:));
plot(best_result_(2,:));
plot(best_result_(3,:));
plot(best_result_(4,:));
set(gca,'FontSize',20);
legend('Double K-means', 'Double NMF', 'Rogers'' Model', 'Coupled NMF');
hold off;

%% Error on Rogers' mothod
err_ = reshape(ERROR_TABLE(3,1,:),[length(set_of_amplitude),1]);
figure;
plot(err_);



