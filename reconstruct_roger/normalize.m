function result = normalize(dataset)
min_ = min(min(dataset));
max_ = max(max(dataset));
result = (dataset-min_)/(max_-min_);
