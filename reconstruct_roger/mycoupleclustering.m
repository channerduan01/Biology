function [Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,HIS] = ...
    MyCoupleClustering(MRNA, PROTEIN, K, J, MAX_ITER, patience, b_verbose)
[T1,N] = size(MRNA);
[T2,~] = size(PROTEIN);
% init parameters
% alternative way to initialize, worse way~
%     Q = rand(K, N);
%     R = rand(J, N, K);
%     for i = 1:N
%         Q(:,i) = Q(:,i)/sum(Q(:,i));
%     end
%     for i = 1:K
%         for j = 1:N
%             R(:,j,i) = R(:,j,i)/sum(R(:,j,i));
%         end
%     end
%     [PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = Maximum(Q,R,MRNA,PROTEIN,K,J,T,N);
THETA = rand(K, J);
for k = 1:K
    THETA(k,:) = THETA(k,:)/sum(THETA(k,:));
end
PI_K = rand(K, 1);
PI_K = PI_K/sum(PI_K);
ii = randperm(N);
AVG_K = MRNA(:,ii(1:K));
VARIANCE_K = ones(K, 1)*mean(var(MRNA,0,2));
ii = randperm(N);
AVG_J = PROTEIN(:,ii(1:J));
VARIANCE_J = ones(J, 1)*mean(var(PROTEIN,0,2));
[Q,R] = Expectation(PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T1,T2,N);

% start optimization
HIS = zeros(MAX_ITER+1,1);
low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T1,T2,N);
if isnan(low_bound), throw(MException('MyCoupleClustering:Nan','Nan value!!!')); end
HIS(1,1) = low_bound;
last_low_bound = low_bound;
for iter = 1:MAX_ITER
    [Q,R] = Expectation(PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T1,T2,N);
    [PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = Maximum(Q,R,MRNA,PROTEIN,K,J,T1,T2,N);
    low_bound = CalcuLowbound(Q,R,PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T1,T2,N);
    if isnan(low_bound), throw(MException('MyCoupleClustering:Nan','Nan value!!!')); end
    HIS(iter+1,1) = low_bound;
    step = low_bound-last_low_bound;
    if b_verbose, fprintf('iter-%d  step: %f, low_bound: %f, dive-R: %f\n', iter, step,...
            low_bound, DivergenceOfR(R,K)); end
    if step > 0 && step < patience
        if b_verbose, fprintf('iter-%d, converged!!!\n', iter); end
        break
    end
    last_low_bound = low_bound;
end
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

function [PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J] = Maximum(Q,R,MRNA,PROTEIN,K,J,T1,T2,N)
PI_K = sum(Q,2)/N;
AVG_K = zeros(T1, K);
for k = 1:K
    sum_ = zeros(T1,1);
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
    VARIANCE_K(k) = sum_/(T1*sum(Q(k,:)));
end
THETA = zeros(K, J);
for k = 1:K
    for j = 1:J
        sum_ = sum(Q(k,:));
        if sum_ ~= 0
            THETA(k, j) = sum(R(j,:,k).*Q(k,:))/sum_;
        else
            THETA(k, j) = 0;
        end
    end
end
AVG_J = zeros(T2, J);
for j = 1:J
    sum_1 = zeros(T2,1);
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
    VARIANCE_J(j) = sum_1/(T2*sum_2);
end
end

function [Q,R] = Expectation(PI_K,AVG_K,VARIANCE_K,THETA,AVG_J,VARIANCE_J,MRNA,PROTEIN,K,J,T1,T2,N)
Q = zeros(K, N);
R = zeros(J, N, K);
TMP_J_SUM_MATRIX = zeros(K, N);
for k = 1:K
    for j = 1:J
        R(j,:,k) = MyMvnpdf(PROTEIN',AVG_J(:,j)',VARIANCE_J(j)*eye(T2))'*THETA(k,j);
    end
    for i = 1:N
        TMP_J_SUM_MATRIX(k,i) = sum(R(:,i,k));
        R(:,i,k) = R(:,i,k)/TMP_J_SUM_MATRIX(k,i);
    end
end
for k = 1:K
    Q(k,:) = MyMvnpdf(MRNA',AVG_K(:,k)',VARIANCE_K(k)*eye(T1))'*PI_K(k).*TMP_J_SUM_MATRIX(k,:);
end
for i = 1:N
    Q(:,i) = Q(:,i)/sum(Q(:,i));
end
end
