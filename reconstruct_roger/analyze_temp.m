%% Roger EM coupled clustering model
%
% init data
% close all
clear
clc

[MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] = GeneDataLoad();
load RESULT;
K = 15;
J = 19;

%%
REPEAT = length(RESULT);
best_idx = 0;
lowest_bound = 0;
for i = 1:REPEAT
    if RESULT{i}.low_bound > lowest_bound
        lowest_bound = RESULT{i}.low_bound;
        best_idx = i;
    end
end

Q = RESULT{best_idx}.Q;
R = RESULT{best_idx}.R;
PI_K = RESULT{best_idx}.PI_K;
AVG_K = RESULT{best_idx}.AVG_K;
VARIANCE_K = RESULT{best_idx}.VARIANCE_K;
THETA = RESULT{best_idx}.THETA;
AVG_J = RESULT{best_idx}.AVG_J;
VARIANCE_J = RESULT{best_idx}.VARIANCE_J;
THETA_reverse = RESULT{best_idx}.THETA_reverse;
R_J = RESULT{best_idx}.R_J;



%% Single Analysis
SELETION_THRESHOLD = 0.3;
[mrna_clusters, protein_clusters] = CalcuClusterExtent(K, J, N, Q, R_J, SELETION_THRESHOLD);
OutputClustersNM(K, J, mrna_clusters, protein_clusters, names);
ANALYSIS_PROTEIN_CLUSTER_IDX = 10;
ClusterAnalysis(ANALYSIS_PROTEIN_CLUSTER_IDX, Q, protein_clusters, MRNA, PROTEIN_ORIGINAL);
 


