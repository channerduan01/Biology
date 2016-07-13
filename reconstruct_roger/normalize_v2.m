function result = normalize_v2(dataset)
result = dataset;
minValues = min(dataset, [], 2);
maxValues = max(dataset, [], 2);
for i = 1:size(dataset,1)
    if minValues(i) < 0
        result(i,:) = result(i,:) - minValues(i);
%         result(i,:) = (result(i,:) - minValues(i))/(maxValues(i)-minValues(i));
    end
end
