function [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,low_bound] = ...
    mycoupleclustering(MRNA, PROTEIN, K, J, MAX_ITER, patience, b_verbose)
    [T,N] = size(MRNA);
    % init expectation
    Q = rand(K, N);
    R = rand(J, N, K);
    for i = 1:N
        Q(:,i) = Q(:,i)/sum(Q(:,i));
    end
    for i = 1:K
        for j = 1:N
            R(:,j,i) = R(:,j,i)/sum(R(:,j,i));
        end
    end
    [PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = Maximum(Q,R,MRNA,PROTEIN,K,J,T,N);
    iter = 0;
    while iter < MAX_ITER
        low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N);
        if iter > 0
            if b_verbose, fprintf('iter-%d  step: %f, dive-R: %f\n', iter, low_bound-last_low_bound, DivergenceOfR(R,K)); end
            if low_bound-last_low_bound < patience
                if b_verbose, fprintf('iter-%d, converged!!!\n', iter); end
                R = mean(R, 3);
                break
            end
        end
        iter = iter+1;
        [Q,R] = Expectation(PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N);
        [PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = Maximum(Q,R,MRNA,PROTEIN,K,J,T,N);
        last_low_bound = low_bound;
    end
    R = mean(R, 3);
end

%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions
%------------------------------------------------------------------------------------------------------------------------
function divergence = DivergenceOfR(R, K)
    meanR = mean(R,3);
    divergence = 0;
    for k = 1:K
        divergence = divergence + sum(sum(abs(R(:,:,k) - meanR)));
    end
end

function low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N)
    item1 = 0;
    item2 = 0;
    item456 = 0;
    LOG_DENSITY_J = zeros(J,N);
    for j = 1:J
        LOG_DENSITY_J(j,:) = log(mvnpdf(PROTEIN',AVG_J(:,j)',VARIANCE_J(j)*eye(T)));
    end
    for i = 1:N
       item1 = item1 + sum(log(PI_K).*Q(:,i));
       item2 = item2 - sum(log(Q(:,i)).*Q(:,i));
       for k = 1:K
           for j = 1:J
               item456 = item456 + Q(k,i)*R(j,i,k)*(log(THETA(k,j)) - log(R(j,i,k)) + LOG_DENSITY_J(j,i));
           end
       end
    end
    item3 = 0;
    for k = 1:K
        item3 = item3 + sum(log(mvnpdf(MRNA',AVG_K(:,k)',VARIANCE_K(k)*eye(T)))'.*Q(k,:));
    end
    low_bound = item1+item2+item3+item456;
end

function [PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = Maximum(Q,R,MRNA,PROTEIN,K,J,T,N)
    PI_K = sum(Q,2)/N;
    AVG_K = zeros(T, K);
    for k = 1:K
        sum_ = zeros(T,1);
        for i = 1:N
            sum_ = sum_ + MRNA(:,i)*Q(k,i);
        end
        AVG_K(:,k) = sum_/sum(Q(k,:));
    end
    VARIANCE_K = zeros(K,1);
    for k = 1:K
        sum_ = 0;
        for i = 1:N
            sum_ = sum_ + sum((MRNA(:,i)-AVG_K(:,k)).^2)*Q(k,i);
        end
        VARIANCE_K(k) = sum_/(T*sum(Q(k,:)));
    end
    THETA = zeros(K, J);
    for k = 1:K
        for j = 1:J
            THETA(k, j) = sum(R(j,:,k).*Q(k,:))/sum(Q(k,:));
        end     
    end
    AVG_J = zeros(T, J);
    for j = 1:J
        sum_1 = zeros(T,1);
        sum_2 = 0;
        for k = 1:K
            for i = 1:N
                sum_1 = sum_1 + PROTEIN(:,i)*R(j,i,k)*Q(k,i);
                sum_2 = sum_2 + R(j,i,k)*Q(k,i);
            end
        end
        AVG_J(:,j) = sum_1/sum_2;
    end
    VARIANCE_J = zeros(J,1);
    for j = 1:J
        sum_1 = 0;
        sum_2 = 0;
        for k = 1:K
            for i = 1:N
                sum_1 = sum_1 + sum((PROTEIN(:,i)-AVG_J(:,j)).^2)*R(j,i,k)*Q(k,i);
                sum_2 = sum_2 + R(j,i,k)*Q(k,i);
            end
        end
        VARIANCE_J(j) = sum_1/(T*sum_2);
    end    
end

function [Q,R] = Expectation(PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N)
    Q = zeros(K, N);
    R = zeros(J, N, K);
    TMP_J_SUM_MATRIX = zeros(K, N);
    for k = 1:K
        for j = 1:J
             R(j,:,k) = mvnpdf(PROTEIN',AVG_J(:,j)',VARIANCE_J(j)*eye(T))'*THETA(k,j);
        end
        for i = 1:N
            TMP_J_SUM_MATRIX(k,i) = sum(R(:,i,k));
            R(:,i,k) = R(:,i,k)/TMP_J_SUM_MATRIX(k,i);
        end
    end
    for k = 1:K
        Q(k,:) = mvnpdf(MRNA',AVG_K(:,k)',VARIANCE_K(k)*eye(T))'*PI_K(k).*TMP_J_SUM_MATRIX(k,:);
    end
    for i = 1:N
        Q(:,i) = Q(:,i)/sum(Q(:,i));
    end
end
