function [RESORT_THETA] = ArrangeTheta(idx_mrna_original,idx_mrna, idx_protein_original, idx_protein, THETA, K, J)

map_mrna = zeros(K,K);
for k1 = 1:K
    for k2 = 1:K
        map_mrna(k1,k2) = length(intersect(find(idx_mrna==k1), find(idx_mrna_original==k2)));
    end
end
map_protein = zeros(J,J);
for j1 = 1:J
    for j2 = 1:J
        map_protein(j1,j2) = length(intersect(find(idx_protein==j1), find(idx_protein_original==j2)));
    end
end
% map_mrna
% map_protein
map_ = map_mrna;
[~,idx_best] = sort(map_(:),1,'descend');
map_tmp_mark = zeros(1, K);
map_res = zeros(1, K);
for i = 1:K*K
    k_row = mod(idx_best(i)-1,K)+1;
    k_column = floor((idx_best(i)-1)/K)+1;
    if (map_tmp_mark(k_row) == 0 && map_res(k_column) == 0)
        map_tmp_mark(k_row) = 1;
        map_res(k_column) = k_row;
    end
    if sum(map_res==0) == 0, break;end
end
% map_res
THETA1 = THETA(map_res,:);

map_ = map_protein;
[~,idx_best] = sort(map_(:),1,'descend');
map_tmp_mark = zeros(1, J);
map_res = zeros(1, J);
for i = 1:J*J
    k_row = mod(idx_best(i)-1,J)+1;
    k_column = floor((idx_best(i)-1)/J)+1;
    if (map_tmp_mark(k_row) == 0 && map_res(k_column) == 0)
        map_tmp_mark(k_row) = 1;
        map_res(k_column) = k_row;
    end
    if sum(map_res==0) == 0, break;end
end
% map_res
RESORT_THETA = THETA1(:,map_res);


