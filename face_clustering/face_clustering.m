%% init
close all
clear
clc

addpath(genpath('/Users/channerduan/Desktop/Final_Project/codes'));

% Load Images from ATT Face Database
faceDatabase = imageSet('/Users/channerduan/Desktop/Final_Project/codes/face_clustering/att_orl_faces', 'recursive');

[images, subject_num, samples_num, size_img, lables] = dataLoad(faceDatabase,true);
N = subject_num*samples_num;
V = normalize(images');
K = 40;


%% K-means clustering
[idx,centres] = kmeans(V',K,'Replicates',1,'Options', statset('Display','final'),'Start','sample');
W = centres';
H = zeros(K,length(idx));
for i = 1:length(idx)
    H(idx(i),i) = 1;
end

displayImages(W, size_img, 3, 'K-means centres', 1:size(H,2));
fprintf('k-means correct rate: %f\n', purity(idx, lables));


%% NMF clustering
% [W,H,~] = mynmf(V,K,'verbose',1,'METHOD','ALS','ALPHA',1,'BETA',1,'MAX_ITER',80,'MIN_ITER',30);
[W,H,~,~,~,~,~] = CoNMF_v2_separate(V, V, K, K ...
        , 'MAX_ITER', 80, 'MIN_ITER', 80, 'VERBOSE', 0, 'METHOD', 'AS' ...
        , 'W_COEF', 1, 'H_COEF', 1, 'T_COEF', 0.5, 'PATIENCE', 0.01 ...
        );
[~, idx] = max(H);
displayImages(W, size_img, 3, 'NMF patterns', 1:size(H,2));
fprintf('NMF correct rate: %f\n', purity(idx, lables));


%% Co-NMF
% create profile matrix
girls_indices = [2,26,29,39];
% girls_indices = [5,6,9,11,12,13,20,21,25,28,31,34,37];
H_info = zeros(2, size(V,2));
for i = 1:size(V,2)
    if any(floor((i-1)/samples_num)+1 == girls_indices)
        H_info(1,i) = 1;
    else
        H_info(2,i) = 1;
    end
end

[W,H,W2_res,H2_res,Theta,~,~] = CoNMF_v4_flow(V, H_info, K, 2 ...
    , 'MAX_ITER', 80, 'MIN_ITER', 20, 'VERBOSE', 0, 'METHOD', 'AS' ...
    , 'W_COEF', 4, 'H_COEF', 4, 'T_COEF', 0.4, 'PATIENCE', 0.01 ...
    );
[~, idx] = max(H);
[~,i] = sort(Theta);
displayImages(W, size_img, 3, 'Co-NMF patterns girls positive', i(:,1));
fprintf('Co-NMF correct rate: %f\n', purity(idx, lables));


%% Reconstruction process
test_img_index = 60;
original_image_size = [112,92];

img = reshape(all_train_img(test_img_index,:), original_image_size);
figure();
colormap(gray);
imagesc(img);
axis('off');
axis('equal');

[comp_weight,comp_idxs] = sort(H(H(:,test_img_index)>0,test_img_index),1,'descend');
number_of_comp = length(comp_idxs);
selected_comp = base_NMF(:,comp_idxs);
reconstruct = zeros(112*92,1);
display_ = zeros(112*92,2*number_of_comp);
for i = 1:number_of_comp
    reconstruct = reconstruct + selected_comp(:,i)*comp_weight(i);
    display_(:,i*2-1) = selected_comp(:,i);
    display_(:,i*2) = normalize(reconstruct);
end

displayImages(display_, [112,92], 3, 'reconstruct process');


%% All subject
subject_num = length(faceDatabase);
display_ = zeros(112*92,subject_num);
for i = 1:subject_num
    display_(:,i) = reshape(read(faceDatabase(i),1), [], 1);
end
displayImages(display_, [112,92], 1, 'subjects', 1:size(display_,1));




