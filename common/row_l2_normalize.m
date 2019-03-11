function out = row_l2_normalize(X) 
out = bsxfun(@times, X, 1./row_l2_norms(X));