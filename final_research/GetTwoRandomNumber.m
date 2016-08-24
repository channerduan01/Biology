%%
% Get two random numbers from 1 to N
%
function [num1, num2] = GetTwoRandomNumber(N)
num1 = floor(rand*N)+1;
num2 = num1;
while(num2 == num1)
    num2 = floor(rand*N)+1;
end
end




