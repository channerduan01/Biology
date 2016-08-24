%%
% Data Generator
%
function [THETA_ORIGINAL] = CreateThetaOriginal(K, J, T, N, MRNA, PROTEIN, AVG_K_, COV_K_, AVG_J_, COV_J_, bShowPattern)
H1_ORIGINAL = zeros(K, N);
H2_ORIGINAL = zeros(J, N);
for k = 1:K
    H1_ORIGINAL(k, :) = MyMvnpdf(MRNA',AVG_K_(:,k)', COV_K_(:,:,k));
    H2_ORIGINAL(k, :) = MyMvnpdf(PROTEIN',AVG_J_(:,k)', COV_J_(:,:,k));
end
for i = 1:N
    H1_ORIGINAL(:, i) = H1_ORIGINAL(:, i)/sum(H1_ORIGINAL(:, i));
    H2_ORIGINAL(:, i) = H2_ORIGINAL(:, i)/sum(H2_ORIGINAL(:, i));
end
THETA_ORIGINAL = CalcuTheta(H1_ORIGINAL, H2_ORIGINAL, K, J, N);
if bShowPattern
    figure();
    imagesc(THETA_ORIGINAL);
    title('Real THETA');
end
end




