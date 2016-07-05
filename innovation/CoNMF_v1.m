function [W1,H1,W2,H2,THETA,HIS,last_cost] = ...
    CoNMF(MRNA, PROTEIN, K, J, wCoef, hCoef, tCoef, MAX_ITER, b_verbose, patience)

[T,N] = size(MRNA);
V1 = MRNA;
V2 = PROTEIN;

% init
W1 = rand(T,K);
H1 = rand(K,N);
% W2 = rand(T,J);
% H2 = rand(J,N);

W2 = zeros(T,J);
H2 = zeros(J,N);


THETA = rand(K,J);
% only for test ===================================
% for i = 1:N
%     H1(:,i) = H1(:,i)/sum(H1(:,i));
%     H2(:,i) = H2(:,i)/sum(H2(:,i));
% end
for k = 1:K
    THETA(k,:) = THETA(k,:)/sum(THETA(k,:));
end

% =================================================


salphaI_K = sqrt(wCoef)*eye(K);
sbetaE_K = sqrt(hCoef)*ones(1,K);
salphaI_J = sqrt(wCoef)*eye(J);
sbetaE_J = sqrt(hCoef)*ones(1,J);
sbetaE_Theta = sqrt(tCoef)*ones(1,K);


zero1N = zeros(1,N);
zero1J = zeros(1,J);
zeroKT = zeros(K,T);
zeroJT = zeros(J,T);

HIS = zeros(MAX_ITER, 5);
last_cost = Inf;
% mode_forward = true;
for i = 1:MAX_ITER
    
    %     if mode_forward
    %     if mod((i-1), 10) < 5
    % first
    %     H1 = updateH([V1;zero1N],[W1;sbetaE_K],H1);
    %     W2 = W1*THETA;
    %     H2 = updateH([V2;zero1N],[W2;sbetaE_J],H2);
    
    % second
    %     W1 = updateW([V1 zeroKT'],W1,[H1 salphaI_K']);
    %     W2 = updateW([V2 zeroJT'],W2,[H2 salphaI_J']);
    %     THETA = updateH([W2;zero1J],[W1;sbetaE_Theta],THETA);
    
    
    %     else
    %             % third
    %             W1 = W2*pinv(THETA);
    %             H1 = updateH([V1;zero1N],[W1;sbetaE_K],H1);
    %             H2 = updateH([V2;zero1N],[W2;sbetaE_J],H2);
    %
    %             % fourth
    %             W1 = updateW([V1 zeroKT'],W1,[H1 salphaI_K']);
    %             W2 = updateW([V2 zeroJT'],W2,[H2 salphaI_J']);
    %             THETA = pinv(updateHTest(W1,W2,pinv(THETA),tCoef));
    %     end
    
    %     %first
    H1 = updateH(V1,W1,H1,hCoef);
    % only for test ===================================
    for ii = 1:N
        H1(:,ii) = H1(:,ii) + (1-sum(H1(:,ii)))/length(H1(:,ii));
    end
    H1(H1<0) = 0;
    % =================================================
    W2 = W1*THETA;
    H2 = updateH(V2,W2,H2,hCoef);
    for ii = 1:N
        H2(:,ii) = H2(:,ii) + (1-sum(H2(:,ii)))/length(H2(:,ii));
    end
    H2(H2<0) = 0;
    
    %     %second
    W1 = updateW(V1,W1,H1,wCoef);
    W2 = updateW(V2,W2,H2,wCoef);
    THETA = updateH(W2,W1,THETA,tCoef);
    for k = 1:K
        THETA(k,:) = THETA(k,:)/sum(THETA(k,:));
    end
    
    %     % third
    %     W1 = W2*pinv(THETA);
    %     H1 = updateH(V1,W1,H1,hCoef);
    %     H2 = updateH(V2,W2,H2,hCoef);
    %     % fourth
    %     W1 = updateW(V1,W1,H1,wCoef);
    %     W2 = updateW(V2,W2,H2,wCoef);
    %     THETA = pinv(updateH(W1,W2,pinv(THETA),tCoef));
    %
    if b_verbose, outputCost(i,V1,W1,H1,V2,W2,H2,THETA);end
    HIS(i,:) = recordCost(V1,W1,H1,V2,W2,H2,THETA);
    promote = last_cost - HIS(i,1);
    last_cost = HIS(i,1);
    %     fprintf('%d,  promote - %f\n', i, promote);
    if promote > 0 && promote < patience
        %        mode_forward = ~mode_forward;
        
        %        last_cost = Inf;
        break;
    end
end
end
%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------
function H = updateH(V,W,H,coef)
H = H .* (W'*V)./(W'*W*H + coef*H + eps);
end
function W = updateW(V,W,H,coef)
W = W .* (V*H')./(W*H*H' + coef*W + eps);
end

% function H = updateH(V,W,H,coef)
% H = pinv(W'*W+coef*eye(size(W,2)))*(W'*V);
% H(H<0) = 0;
% end
%
% function W = updateW(V,W,H,coef)
% W = (pinv(H*H'+coef*eye(size(H,1)))*(H*V'))';
% W(W<0) = 0;
% end

% function H = updateHTest(V,W,H,coef)
%     H = H .* (W'*V)./(W'*W*H + coef*H + eps);
% end

% function H = updateH(V,W,H)
% [H,~,~] = nnlsm_blockpivot(W,V,0,H);
% end
% function W = updateW(V,W,H)
% [W,~,~] = nnlsm_blockpivot(H',V',0,W');
% W = W';
% end

% function H = updateH(V,W,H)
% [H,~,~] = nnlsm_activeset(W,V,1,0,H);
% end
% function W = updateW(V,W,H)
% [W,~,~] = nnlsm_activeset(H',V',1,0,W');
% W = W';
% end

function c = cost(V,W,H)
c = sqrt(sum(sum((V-W*H).^2)));
end

function outputCost(i,V1,W1,H1,V2,W2,H2,THETA)
record = recordCost(V1,W1,H1,V2,W2,H2,THETA);

fprintf('iter-%d  Whole-cost:%f, V1-cost: %f, V2-cost: %f, CO-cost: %f - %f\n', ...
    i, record(1), record(2), record(3), record(4), record(5));
end

function record = recordCost(V1,W1,H1,V2,W2,H2,THETA)
a = cost(V1,W1,H1);
b = cost(V2,W2,H2);
c = cost(W2,W1,THETA);
d = cost(W1,W2,pinv(THETA));
record = [a+b+c+d a b c d];
end





