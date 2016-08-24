%%
% Get two random numbers from 1 to N
%
function [resList] = GetRandomList(basicList, number)
resList = zeros(1, number);
size_ = length(basicList);
for i = 1:number
    idx_ = floor(rand*size_)+1;
    resList(i) = basicList(idx_);
    basicList(idx_) = basicList(size_);
    size_ = size_-1;
end
end




