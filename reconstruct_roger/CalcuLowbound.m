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
               item456 = item456 + Q(k,i)*R(j,i,k)*(log(THETA(k,j)) - log(R(j,i,k))+ LOG_DENSITY_J(j,i));
           end
       end
    end
    item3 = 0;
    for k = 1:K
        item3 = item3 + sum(log(mvnpdf(MRNA',AVG_K(:,k)',VARIANCE_K(k)*eye(T)))'.*Q(k,:));
    end
    low_bound = item1+item2+item3+item456;
end