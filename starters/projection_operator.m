% The component for Project gradient descent
% inspired by Hoyer's idea in Non-negative Matrix Factorization with Sparseness Constraints

function [s, iter] = projection_operator(x,L1,L2)

iter = 0;
dim = length(x);
s = x + (L1-sum(x))/dim;
z = zeros(dim,1);
while true
    iter = iter + 1;
    m = zeros(dim,1);
    m(~z) = L1/(dim-sum(z));
    if sum(s==m) == dim
        break;
    end
%     eq = [num2str(sum((s-m).^2)),'*a^2+',num2str(sum(m.*(s-m))),'*a+',num2str(sum(m.^2)),'-',num2str(L2^2),'=0'];
%     quadr = solve(eq,'a');
    quadr = roots([sum((s-m).^2) sum(m.*(s-m)) sum(m.^2)-L2^2]);
    a = quadr(quadr>0);
    if length(a) ~= 1
        error('Crash! Unexpected quadratic results for a! (all negative)');
    end
    s = m + a*(s-m);
    if sum(s<0) == 0
        break;
    end
    z = z + s<0;
    s(s<0) = 0;
    s(~z) = s(~z) - (sum(s)-L1)/(dim-sum(z));
end

end