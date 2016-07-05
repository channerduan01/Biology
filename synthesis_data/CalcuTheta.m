function [THETA_ORIGINAL] = CalcuTheta(H1, H2, K, J, N)

PI_K_ = sum(H1,2)/N;
THETA_ORIGINAL = zeros(K, J);
for i = 1:N
    THETA_ORIGINAL = THETA_ORIGINAL + H1(:, i)*H2(:, i)';
end
THETA_ORIGINAL = THETA_ORIGINAL ./ sum(sum(THETA_ORIGINAL));
for j = 1:J
    THETA_ORIGINAL(:,j) = THETA_ORIGINAL(:,j)./PI_K_;
end


