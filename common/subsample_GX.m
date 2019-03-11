function [A, l1, l2, X_1, X_2] = subsample_GX(alpha, A, l1, l2, X_1, X_2)

N = size(A);
idx1 = randsample(N(1), round(alpha*N(1)));
idx2 = randsample(N(2), round(alpha*N(2)));

A = A(idx1,idx2);

l1 = l1(idx1);
l2 = l2(idx2);

if ~isempty(X_1)
    X_1 = X_1(idx1,:);
end

if ~isempty(X_2)
    X_2 = X_2(idx2,:);
end