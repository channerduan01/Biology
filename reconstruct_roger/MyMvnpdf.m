%% tricky function, for scaling problem
function [res_] = MyMvnpdf(data, mean, cov)
if sum(sum(isnan(cov))) > 0
    res_ = zeros(size(data,1),1);
    return;
end

activate_ = diag(cov) ~= 0;
dimension = sum(activate_);
switch dimension
    case 0
        res_ = zeros(size(data,1),1);
    case 1
        res_ = normpdf(data(:,activate_), mean(:,activate_), cov(activate_,activate_));
    otherwise
        res_ = mvnpdf(data(:,activate_), mean(:,activate_), cov(activate_,activate_));
end