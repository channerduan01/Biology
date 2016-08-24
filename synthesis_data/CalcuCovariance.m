function [C_] = CalcuCovariance(X)
[N, M] = size(X);
C_ = zeros(M, M);
mean_row = mean(X);
for i = 1:N
    C_ = C_ + (X(i,:)-mean_row)'*(X(i,:)-mean_row);
end
C_ = C_./N;

