function [A,noise_A] = addnoise_GX(lambda, A)

N = size(A);

p = lambda*sum(N)/(2*prod(N));
noise_A = sparse(rand(N(1),N(2)) < p  );
A = A + noise_A;
