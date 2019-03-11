function [ A ] = genSBM3( B, y, z )
%GENDCBLKMOD Summary of this function goes here
%   Detailed explanation goes here
%   we will assume y and z are in matrix form

ny = sum(y);
nz = sum(z);
y = label_mat2vec(y);
z = label_mat2vec(z);

%[K, L] = size(B);
% A = zeros(sum(ny), sum(nz));
% for k = 1:K
%     for ell = 1:L
%         idx_y = y == k;
%         idx_z = z == ell;
%         A( idx_y , idx_z ) = binornd(1, B(k,ell), [ny(k), nz(ell)]); 
%     end
% end
% A = sparse(A);



%function As = KblockGraph(csizes, Pmat)
% K blocks generation
  % csizes: vector of block sizes
  % Pmat: density within/between blocks

  indexrr = [];
  indexcc = [];
  %cumcsizes = cumsum(csizes);
  cumny = cumsum(ny);
  cumnz = cumsum(nz);
  K = length(ny);
  L = length(nz);
  n = cumny(K);
  m = cumnz(L);
  for k=1:K,
    for ell=1:L,
        
      [tmprr, tmpcc] = offblock(ny(k), nz(ell), B(k,ell));
          
      if (k > 1)
        tmprr = tmprr + cumny(k - 1);
      end
      if (ell > 1)
        tmpcc = tmpcc + cumnz(ell - 1);
      end
      indexrr = [indexrr, tmprr];
      indexcc = [indexcc, tmpcc];
      
    end
  end
  % adjust the order according to "y" and "z"
  [~, Iy] = sort(y);
  [~, Iz] = sort(z);
  ii = Iy(indexrr);
  jj = Iz(indexcc);
  A = sparse(ii, jj, 1, n, m);
end


function [rr, cc] = offblock(nr, nc, p) 
  rr=[]; 
  cc=[];
  if (p > 0) 
      nsize = nr * nc;
      if (nsize < 1e7)
        msample = binornd(nsize, p, 1, 1);
      elseif (nsize*p > 1e6),
        msample = round(normrnd(nsize*p, sqrt(nsize*p*(1-p)),1,1)); 
        msample = max(0, msample);
      else
        msample = poissrnd(nsize*p, 1, 1);
      end
      if (msample > 0)
	    index = mysamplefun(nsize, msample);
        rr = ceil(index/nc);
        cc = index - (rr-1) * nc;
      end
  end
end

