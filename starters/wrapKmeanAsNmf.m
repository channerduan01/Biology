function [W,H] = wrapKmeanAsNmf(A,k)
[idx,centres] = kmeans(A',k,'Replicates',1,'Start','sample');
W = centres';
H = zeros(k,length(idx));
for i = 1:length(idx)
    H(idx(i),i) = 1;
end
end