function s = compDevFromMean(A,y,z,yP,zPt)
Af = struct;
Af.size = size(A);
Af.times = @(x,varargin) A*x - yP * (z'*x) ;
Af.trans = @(x,varargin) A'*x - zPt * (y'*x);
Af.param = [];
opts.tol = 1e-8;
[~,s,~] = lmsvd(Af,1,[]);