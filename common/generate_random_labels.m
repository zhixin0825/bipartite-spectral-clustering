function out = generate_random_labels(n,temp)

if length(temp) == 1
    K = temp;
    prior = ones(1,K)/K;
else
    prior = temp(:)';
    K = length(temp);
    prior = prior/sum(prior);
end

%%
out = mnrnd(1,prior,n);
% X = repmat(rand(n,1),1,K);
% thresh = [zeros(n,1) repmat(cumsum(prior),n,1)];
% temp = X < thresh;
% [X temp]
% diff(temp,1,2)