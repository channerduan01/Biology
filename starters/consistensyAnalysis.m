% Using connectivity matrix to show the consistency of an NMF clustering
% algorithm
function [pk,COST_,C_] = consistensyAnalysis(A,range,repeatTime,fun_nmf)
cost = @(A,W,H) sqrt(sum(sum((A-W*H).^2)));
size = length(A);
number_k = length(range);
C_ = zeros(number_k,size,size);
C_TMP = zeros(repeatTime,size^2);
COST_ = zeros(number_k,repeatTime);
pk = zeros(number_k,1);
for i = 1:number_k
    k = range(i);
    for j = 1:repeatTime
        [W,H] = fun_nmf(A,k);
        [~,idx] = max(H);
        C_TMP(j,:) = reshape(getConnectivityMatrix(idx),1,size^2);
        COST_(i,j) = cost(A,W,H);
    end
    C_(i,:,:) = reshape(mean(C_TMP),size,size);
    pk(i) = 4/size^2*sum(sum((C_(i,:,:)-0.5).^2));
end
end

%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------
function C = getConnectivityMatrix(idx)
size = length(idx);
C = zeros(size,size);
for i = 1:size
    C(i,:) = idx==idx(i);
end
end