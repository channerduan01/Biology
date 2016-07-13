function result = normalize(dataset)
% global way
min_ = min(min(dataset));
max_ = max(max(dataset));
result = dataset-min_;
% result = (dataset-min_)/(max_-min_);
