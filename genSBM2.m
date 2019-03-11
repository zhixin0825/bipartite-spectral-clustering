function [ A ] = genSBM2( B, y, z )
%GENDCBLKMOD Summary of this function goes here
%   Detailed explanation goes here
%   we will assume y and z are in matrix form

ny = sum(y);
nz = sum(z);
y = label_mat2vec(y);
z = label_mat2vec(z);

[K, L] = size(B);
A = zeros(sum(ny), sum(nz));
for k = 1:K
    for ell = 1:L
        idx_y = y == k;
        idx_z = z == ell;
        A( idx_y , idx_z ) = binornd(1, B(k,ell), [ny(k), nz(ell)]); 
    end
end
A = sparse(A);

