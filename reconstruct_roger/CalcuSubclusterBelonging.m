function [R_J, Q_J] = CalcuSubclusterBelonging( ...
    MRNA, AVG_K, VARIANCE_K, PROTEIN, AVG_J, VARIANCE_J, PI_K, THETA_reverse, R, K, J, N, T1, T2)
R_J = zeros(J, N);
for j = 1:J
    for i = 1:N
        tmp = R(j,i,:);
        R_J(j,i) = sum(tmp(:).*PI_K);
    end
end
p_tmp_mrna = zeros(N, K);
for k = 1:K
    p_tmp_mrna(:, k) = MyMvnpdf(MRNA',AVG_K(:,k)',VARIANCE_K(k)*eye(T1));
end
p_tmp_protein = zeros(N, J);
for j = 1:J
    p_tmp_protein(:, j) = MyMvnpdf(PROTEIN',AVG_J(:,j)',VARIANCE_J(j)*eye(T2));
end
PI_J = sum(R_J,2)/N;
Q_J = zeros(J, N);
for j = 1:J
    for i = 1:N
        for k = 1:K
            Q_J(j,i) = Q_J(j,i) + PI_J(j)*THETA_reverse(j,k)*p_tmp_mrna(i,k)*p_tmp_protein(i,j);
        end
    end
end
for i = 1:N
    Q_J(:,i) = Q_J(:,i)/sum(Q_J(:,i));
end