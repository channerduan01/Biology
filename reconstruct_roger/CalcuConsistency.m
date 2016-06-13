% Input IDX_MATRIX: 
%   rows number is repeat times,
%   columns number is items number of basic data
%
function [pk,C_] = CalcuConsistency(IDX_MATRIX)
    [repeatTime, N] = size(IDX_MATRIX);
    C_TMP = zeros(repeatTime,N^2);
    for i = 1:repeatTime
        C_TMP(i,:) = reshape(getConnectivityMatrix(IDX_MATRIX(i,:)),1,N^2);
    end
    C_ = reshape(mean(C_TMP),N,N);
    pk = 4/N^2*sum(sum((C_-0.5).^2));
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