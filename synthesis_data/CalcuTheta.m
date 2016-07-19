function [THETA] = CalcuTheta(H1, H2, K, J, N)
H1 = normalizeColumn(H1);
H2 = normalizeColumn(H2);
PI_K_ = sum(H1,2)/N;
THETA = H1*H2';
THETA = THETA ./ sum(sum(THETA));
for j = 1:J
    THETA(:,j) = THETA(:,j)./PI_K_;
end


function A = normalizeColumn(A)
for i = 1:size(A,2)
    sum_ = sum(A(:,i));
    % if the sum is valid and not close to 1
    if sum_ ~=0 && ( sum_ < 0.999 || sum_ > 1.001 )
        A(:,i) = A(:,i)/sum(A(:,i));
    end
end