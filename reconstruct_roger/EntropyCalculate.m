function [THETA_reverse, entropy_j_k, entropy_k_j] = BasicCoupleClusteringAnalysis(K, J, PI_K, THETA)

    
    
    THETA_reverse = zeros(fliplr(size(THETA)));
    for j = 1:J
        P_j = sum(THETA(:,j).*PI_K);
        for k = 1:K
            THETA_reverse(j,k) = THETA(k,j)*PI_K(k)/P_j;
        end
    end
    entropy_j_k = 0;
    for k = 1:K
        for j = 1:J
            entropy_j_k = entropy_j_k + THETA(k,j)*log2(THETA(k,j));
        end
    end
    entropy_j_k = -1/K*entropy_j_k;
    entropy_k_j = 0;
    for j = 1:J
        for k = 1:K
            entropy_k_j = entropy_k_j + THETA_reverse(j,k)*log2(THETA_reverse(j,k));
        end
    end
    entropy_k_j = -1/J*entropy_k_j;
end
