function [THETA] = CalcuTheta(H1, H2, K, J, N)
PI_K_ = sum(H1,2)/N;
THETA = H1*H2';
THETA = THETA ./ sum(sum(THETA));
for j = 1:J
    THETA(:,j) = THETA(:,j)./PI_K_;
end
