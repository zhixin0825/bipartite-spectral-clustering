function [ y ] = truelabel( n )
%TRUELABEL Summary of this function goes here
%   Detailed explanation goes here
y = [];
K = size(n,2);
for k = 1:K
    y = vertcat(y,k*ones(n(k),1));
end
y = ytoPi(y,K);

end

