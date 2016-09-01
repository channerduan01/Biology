%%
% Methods comparison

close all;
clear;
clc;
addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

T = 10;
N = 1000;
K = 5;
J = 5;

[MRNA_ORIGINAL, PROTEIN_ORIGINAL, H1_ORIGINAL__, H2_ORIGINAL__, AVG_K_, COV_K_, AVG_J_, COV_J_] = GenerateBasisData(K, J, T, N);
IDEAL_THETA_ORIGINAL = CalcuTheta(H1_ORIGINAL__, H2_ORIGINAL__, K, J, N);
THETA_ORIGINAL = CreateThetaOriginal(K, J, T, N, MRNA_ORIGINAL, PROTEIN_ORIGINAL, AVG_K_, COV_K_, AVG_J_, COV_J_, true);

%%
MRNA_ORIGINAL = MRNA;
PROTEIN_ORIGINAL = PROTEIN;
H1_ORIGINAL__ = H1_ORIGINAL;
H2_ORIGINAL__ = H2_ORIGINAL;

%%
% create specific noise
NOISE_1 = randn(size(MRNA_ORIGINAL));
NOISE_2 = randn(size(PROTEIN_ORIGINAL));

WHOLE_REPEAT_NUM = 1;
REPEAT_NUM = 5;
set_of_amplitude = 0:0.2:1;
% set_of_amplitude = 0:0.02:0.4;

% set_of_amplitude = [0, 1]

HUGE_TALBE = zeros(REPEAT_NUM, 4*5, length(set_of_amplitude), WHOLE_REPEAT_NUM);
BEST_HUGE_TALBE = zeros(4, 3, length(set_of_amplitude), WHOLE_REPEAT_NUM);
ERROR_TABLE = zeros(4, 1, length(set_of_amplitude), WHOLE_REPEAT_NUM);



%%

for tt = 1:WHOLE_REPEAT_NUM
ii = randperm(N);
for i = 1:length(set_of_amplitude)
    fprintf('\n\n\n>>>>> deal with amplitude: %f\n', set_of_amplitude(i));
    MRNA = normalize(MRNA_ORIGINAL);
    PROTEIN = normalize(PROTEIN_ORIGINAL);
    H1_ORIGINAL = H1_ORIGINAL__;
    H2_ORIGINAL = H2_ORIGINAL__;

    disorder_length = ceil(length(ii)*set_of_amplitude(i));
%     disorder_length = N;    % for testing!
    
    % ------------ I have a good idea to prove the property of them
    % activate this block or not!!!
%     MRNA(:,ii(1:disorder_length)) = MRNA(:,1:disorder_length);
%     H1_ORIGINAL(:,ii(1:disorder_length)) = H1_ORIGINAL(:,1:disorder_length);    
    %-------------
    PROTEIN(:,ii(1:disorder_length)) = PROTEIN(:,1:disorder_length);
    H2_ORIGINAL(:,ii(1:disorder_length)) = H2_ORIGINAL(:,1:disorder_length);

    THETA_ORIGINAL = CreateThetaOriginal(K, J, T, N, MRNA, PROTEIN, AVG_K_, COV_K_, AVG_J_, COV_J_, true);
    % add noise ============================== !!!!!!!!!!!!!!!!!!!!!
    amptitude_noise = 0.5;
%     amptitude_noise = set_of_amplitude(i);
%     MRNA(MRNA>0) = MRNA(MRNA>0) + amptitude_noise*NOISE_1(MRNA>0);
%     PROTEIN(PROTEIN>0) = PROTEIN(PROTEIN>0) + amptitude_noise*NOISE_2(PROTEIN>0);
    % dramatical change!!!
%     MRNA = MRNA + amptitude_noise*NOISE_1;
%     PROTEIN = PROTEIN + amptitude_noise*NOISE_2;
    MRNA(MRNA_ORIGINAL>0) = MRNA(MRNA_ORIGINAL>0) + amptitude_noise*NOISE_1(MRNA_ORIGINAL>0);
    PROTEIN(PROTEIN_ORIGINAL>0) = PROTEIN(PROTEIN_ORIGINAL>0) + amptitude_noise*NOISE_2(PROTEIN_ORIGINAL>0);
    % ========================================
    MRNA = normalize(MRNA);
    PROTEIN = normalize(PROTEIN);
    % run algorithms
    [HUGE_TALBE(:,:,i,tt), BEST_HUGE_TALBE(:,:,i,tt), ERROR_TABLE(:,:,i,tt)] ...
        = RunAllAlgorithm(REPEAT_NUM, MRNA, PROTEIN, THETA_ORIGINAL, H1_ORIGINAL, H2_ORIGINAL, K, J, T, N, false ...
        , 5);
end

end

%% RESULT
result_ = reshape(mean(mean(HUGE_TALBE,4)), [20,length(set_of_amplitude)]);
figure;
axis_x = set_of_amplitude;
%----
% for i = 1:length(axis_x)
%     axis_x(i) = snr(MRNA_ORIGINAL, axis_x(i)*NOISE_1);
% end
title('Effect of Noise on accuracy', 'FontSize', 20)
%----
hold on;
plot(axis_x, result_(1,:)*100);
plot(axis_x, result_(6,:)*100);
plot(axis_x, result_(11,:)*100);
plot(axis_x, result_(16,:)*100);
set(gca,'FontSize',20);
legend('Double K-means', 'Double NMF', 'Rogers'' Model', 'Coupled NMF');
hold off;
%% RESULT
result_ = reshape(mean(mean(HUGE_TALBE,4)), [20,length(set_of_amplitude)]);
figure;
axis_x = set_of_amplitude;
%----
% for i = 1:length(axis_x)
%     axis_x(i) = snr(PROTEIN_ORIGINAL, axis_x(i)*NOISE_2);
% end
title('Effect of Noise on accuracy', 'FontSize', 20)
%----
hold on;
plot(axis_x, result_(2,:)*100);
plot(axis_x, result_(7,:)*100);
plot(axis_x, result_(12,:)*100);
plot(axis_x, result_(17,:)*100);
% title('Effect of Disorder on accuracy', 'FontSize', 20)
set(gca,'FontSize',20);
legend('Double K-means', 'Double NMF', 'Rogers'' Model', 'Coupled NMF');
hold off;
%% RESULT
result_ = reshape(mean(mean(HUGE_TALBE,4)), [20,length(set_of_amplitude)]);

figure;
hold on;
plot(set_of_amplitude, result_(3,:)*100);
plot(set_of_amplitude, result_(8,:)*100);
plot(set_of_amplitude, result_(13,:)*100);
plot(set_of_amplitude, result_(18,:)*100);
title('Effect of Noise on accuracy', 'FontSize', 20)
% title('Effect of Disorder on accuracy', 'FontSize', 20)
set(gca,'FontSize',20);
legend('Double K-means', 'Double NMF', 'Rogers'' Model', 'Coupled NMF');
hold off;

%% Best RESULT

best_ = mean(BEST_HUGE_TALBE, 4);
best_result_ = reshape(best_(:,1,:), [4,length(set_of_amplitude)]);

figure;
hold on;
plot(set_of_amplitude, best_result_(1,:));
plot(set_of_amplitude, best_result_(2,:));
plot(set_of_amplitude, best_result_(3,:));
plot(set_of_amplitude, best_result_(4,:));
set(gca,'FontSize',20);
legend('Double K-means', 'Double NMF', 'Rogers'' Model', 'Coupled NMF');
hold off;

%% Error on Rogers' mothod
error_ = mean(ERROR_TABLE, 4);
err_ = reshape(error_(3,1,:),[length(set_of_amplitude),1]);
figure;
plot(err_);



