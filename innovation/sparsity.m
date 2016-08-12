%% sparsity meaturement
function result = sparsity(v)

[m, n] = size(v);
result = zeros(1,n);
for i = 1:n
    result(i) = sparse(v(:,i),m);
end

result = result(~isnan(result));

end

function value = sparse(v, m)
    value = (sqrt(m)-norm(v,1)/norm(v,'fro'))/(sqrt(m)-1);
end

