%% clustering purity metrics
function purity_value = purity(idx, lables)
N = length(idx);
values = unique(idx);
count_ = 0;
for i = 1:length(values)
    select_ = idx == values(i);
    count_ = count_ + max(histc(lables(select_),unique(lables(select_))));
end
purity_value = count_/N;

