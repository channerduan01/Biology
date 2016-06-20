function [W1,H1,W2,H2,THETA,HIS] = ...
    CoNMF(MRNA, PROTEIN, K, J, hCoef, wCoef, tCoef, MAX_ITER, b_verbose)

[T,N] = size(MRNA);
V1 = MRNA;
V2 = PROTEIN;

% init
W1 = rand(T,K);
H1 = rand(K,N);
W2 = rand(T,J);
H2 = rand(J,N);
THETA = rand(K,J);

HIS = zeros(MAX_ITER, 4);

for i = 1:MAX_ITER
    if b_verbose, outputCost(i,V1,W1,H1,V2,W2,H2,THETA);end
    
    %     if rand < 0.5
%     if mod((i-1), 20) < 10
        % first
        H1 = updateH(V1,W1,H1,hCoef);
%         W2 = W1*THETA;
%         H2 = updateH(V2,W2,H2,hCoef);
        % second
        W1 = updateW(V1,W1,H1,wCoef);
%         W2 = updateW(V2,W2,H2,wCoef);
%         THETA = updateH(W2,W1,THETA,tCoef);
%     else
%         % third
%         W1 = W2*pinv(THETA);
%         H1 = updateH(V1,W1,H1,hCoef);
%         H2 = updateH(V2,W2,H2,hCoef);
%         % fourth
%         W1 = updateW(V1,W1,H1,wCoef);
%         W2 = updateW(V2,W2,H2,wCoef);
%         THETA = pinv(updateH(W1,W2,pinv(THETA),tCoef));
%     end
%     
    HIS(i,:) = recordCost(V1,W1,H1,V2,W2,H2,THETA);
end
end
%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------
% function H = updateH(V,W,H,coef)
%     H = H .* (W'*V)./(W'*W*H + coef*H + eps);
% end
% function W = updateW(V,W,H,coef)
%     W = W .* (V*H')./(W*H*H' + coef*W + eps);
% end

% function H = updateH(V,W,H,coef)
% H = pinv(W'*W+coef*eye(size(W,2)))*(W'*V);
% H(H<0) = 0;
% end
% 
% function W = updateW(V,W,H,coef)
% W = (pinv(H*H'+coef*eye(size(H,1)))*(H*V'))';
% W(W<0) = 0;
% end

function H = updateH(V,W,H,coef)
[H,~,~] = nnlsm_activeset(W,V,1,0,H);
end
function W = updateW(V,W,H,coef)
[W,~,~] = nnlsm_activeset(H',V',1,0,W');
W = W';
end

function c = cost(V,W,H)
c = sqrt(sum(sum((V-W*H).^2)));
end

function outputCost(i,V1,W1,H1,V2,W2,H2,THETA)
fprintf('iter-%d   V1-cost: %f, V2-cost: %f, CO-cost: %f - %f\n', ...
    i, cost(V1,W1,H1), cost(V2,W2,H2), cost(W2,W1,THETA), cost(W1,W2,pinv(THETA)));
end

function c = recordCost(V1,W1,H1,V2,W2,H2,THETA)
c = [cost(V1,W1,H1), cost(V2,W2,H2), cost(W2,W1,THETA), cost(W1,W2,pinv(THETA))];
end





