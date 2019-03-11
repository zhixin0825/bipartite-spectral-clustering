function out = row_normalize_ell2(X)

out = bsxfun(@rdivide, X, sum(X.^2,2).^0.5);