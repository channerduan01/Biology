function low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T,N)
switch_on_protect = false;   % tricky, for scaling problem
item1 = 0;
item2 = 0;
item456 = 0;
LOG_DENSITY_J = zeros(J,N);
for j = 1:J
    if switch_on_protect
        LOG_DENSITY_J(j,:) = MyMvnpdf(PROTEIN',AVG_J(:,j)',VARIANCE_J(j)*eye(T));
    else
        LOG_DENSITY_J(j,:) = log(MyMvnpdf(PROTEIN',AVG_J(:,j)',VARIANCE_J(j)*eye(T)));
    end
end
for i = 1:N
    if switch_on_protect
        item1 = item1 + sum(PI_K.*Q(:,i));
        item2 = item2 - sum(Q(:,i).*Q(:,i));
    else
        item1 = item1 + sum(log(PI_K).*Q(:,i));
        item2 = item2 - sum(log(Q(:,i)).*Q(:,i));
    end
    
    
    for k = 1:K
        for j = 1:J
            if switch_on_protect
                item456 = item456 + Q(k,i)*R(j,i,k)*( THETA(k,j)-R(j,i,k)+LOG_DENSITY_J(j,i));
            else
                item456 = item456 + Q(k,i)*R(j,i,k)*(log(THETA(k,j)) - log(R(j,i,k)) + LOG_DENSITY_J(j,i));
            end
        end
    end
end
item3 = 0;
for k = 1:K
    if switch_on_protect
        item3 = item3 + sum(MyMvnpdf(MRNA',AVG_K(:,k)',VARIANCE_K(k)*eye(T))'.*Q(k,:));
    else
        item3 = item3 + sum(log(MyMvnpdf(MRNA',AVG_K(:,k)',VARIANCE_K(k)*eye(T)))'.*Q(k,:));
    end
end

low_bound = item1+item2+item3+item456;
if isnan(low_bound)
    fprintf('!!! test: %d %d %d %d\n', item1, item2, item3, item456);
end
end





