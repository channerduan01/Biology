function A = normalizeColumn(A)
for i = 1:size(A,2)
    sum_ = sum(A(:,i));
    if sum_ ~= 0 && (sum_ < 0.999 || sum_ > 1.001)
        A(:,i) = A(:,i)/sum(A(:,i));
    end
end
end