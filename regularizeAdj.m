function A = regularizeAdj(A, varargin)

parser = inputParser;
addOptional(parser,'dmax1',nan)  
addOptional(parser,'dmax2',nan)  
addOptional(parser,'tau', 1)   % regularization parameter
addOptional(parser,'ell1', true)   % regularization parameter
parse(parser, varargin{:});

ell1 = parser.Results.ell1; 
tau = parser.Results.tau;
dmax1 = parser.Results.dmax1;
dmax2 = parser.Results.dmax2;

if tau == Inf 
   warning('no regularization (c=Inf).')
   return
end

n = size(A);
    
%row regularization
d1 = sum(A,2);
if isnan(dmax1)
    db1 = mean(d1);
    d_sorted = sort(d1,'descend');
    alpha1 = floor(n(1)/db1);
    thresh1 = tau * d_sorted(alpha1);
else
    thresh1 = tau * dmax1;
end

%column regulariztion
d2 = sum(A,1);
if isnan(dmax2)
    db2 = mean(d2);
    d_sorted = sort(d2,'descend');
    alpha2 = floor(n(2)/db2);
    thresh2 = tau * d_sorted(alpha2);
else
    thresh2 = tau * dmax2;
end

idx1 = d1 > thresh1;
idx2 = d2 > thresh2;

if sum(idx1) == 0 
    warning('no regularization.')
else
    if ell1
        A(idx1,:) = diag(thresh1./d1(idx1)) * A(idx1,:);
    else
        A(idx1,:) = diag(sqrt(thresh1./d1(idx1))) * A(idx1,:);
    end
end
if sum(idx2) == 0 
    warning('no regularization.')
else
    if ell1
        A(:,idx2) = A(:,idx2) * diag(thresh2./d2(idx2)) ;
    else
        A(:,idx2) = A(:,idx2) * diag(sqrt(thresh2./d2(idx2))) ;
    end
end

   
end