function [label_1 ,label_2, Z_2]= biSpecClust3(A,KL,varargin)

% KL: could be 1x1 or 1x2
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'matched',true)
addOptional(parser,'normalize',false) % whether to noramlize rows of the singular vector matrix
addOptional(parser,'perturb',false)
addOptional(parser,'tau',0.1) % default perturbation parameter
addOptional(parser,'pert_geom',true) % geometric versus arithmatic scaling
addOptional(parser,'kmeans_rep',10) 
addOptional(parser,'type','U')  % could be 'U' or 'US', defualt is U
addOptional(parser,'rmax', min(KL))  % maximum number of singular values kept
addOptional(parser,'idx',1:min(KL))
%addOptional(parser,'LAPLACIAN', true) 

parse(parser, varargin{:});
matched = parser.Results.matched;
normalize = parser.Results.normalize;
perturb = parser.Results.perturb;
tau = parser.Results.tau;
pert_geom = parser.Results.pert_geom;
kmeans_rep = parser.Results.kmeans_rep;
type = parser.Results.type;
K = parser.Results.rmax;
idx = parser.Results.idx;

N = size(A);

D_1 = sum(A,2);
D_2 = sum(A,1);

if pert_geom % scaling of the perturbation
    taus = tau/sqrt(prod(N)); %geometric mean 
else
    taus = tau/(sum(N)/2); %arithmatic mean 
end

if ~perturb
    g1 = safe_diag_pwr(D_1);
    g2 = safe_diag_pwr(D_2);
    A_n = diag(g1)*A*diag(g2);
    [U,~,V] = svds(A_n,K);
else
    g1 = safe_diag_pwr(D_1(:) + taus*N(2)*ones(N(1),1));
    g2 = safe_diag_pwr(D_2(:) + taus*N(1)*ones(N(2),1));
    
    %[U,~,V] = svds(@(x,tflag) Afun(x,tflag,A,g1,g2,taus), N, K+1);
    Af = struct;
    Af.size = N;
    Af.times = @(x,varargin) bsxfun(@times,g1,A*(bsxfun(@times,g2,x)))+ taus*g1*(g2'*x); %Afun(x,'notransp',A,g1,g2,taus);
    Af.trans = @(x,varargin) bsxfun(@times,g2,A'*bsxfun(@times,g1,x)) + taus*g2*(g1'*x); %Afun(x,'transp',A,g1,g2,taus);
    Af.param = [];
    opts.tol = 1e-8;
    [U,~,V,] = lmsvd(Af,K,opts);
    %[U,~,V] = lansvd(@(x) Atfun(x,A*1.,g1,g2,taus), @(x) Atfun(x,A'*1.,g2,g1,taus), N(1), N(2), K+1);
end

U_2 = U(:,idx);
V_2 = V(:,idx);
    
    %row_l2_norms = @(X) sqrt(sum(X.^2,2));
    %row_l2_normalize = @(X) bsxfun(@times, X, 1./row_l2_norms(X));
if normalize
    U_2=row_l2_normalize(U_2);
    V_2=row_l2_normalize(V_2);
end


if matched % matched clustering
    %[e, C] = kmeans(Z_2,K,'Replicates',50, 'Start','sample');
    Z_2 = [U_2;V_2];
    [e, C] = kmeans(Z_2,K,'Replicates', kmeans_rep);
    label_1 = e(1:N(1));
    label_2 = e(N(1) + (1:N(2)));
else % unmatched clustering
    switch lower(type)
        case 'u' % default
            label_1 = kmeans(U, KL(1), 'replicates', kmeans_rep, 'onlinephase','off');%,'Options',kmopts);
            label_2 = kmeans(V, KL(2), 'replicates', kmeans_rep, 'onlinephase','off');%,'Options',kmopts);
        case 'us'
            label_1 = kmeans(U*S, KL(1), 'replicates', kmeans_rep,'onlinephase','off');%,'Options',kmopts);
            label_2 = kmeans(V*S, KL(2), 'replicates', kmeans_rep,'onlinephase','off');%,'Options',kmopts);
    end
end
%figure(1), clf, scatter3(Z_2(:,1),Z_2(:,2),Z_2(:,3),'.'), hold on, scatter3(C(:,1), C(:,2), C(:,3),'ro')
        
end

function g = safe_diag_pwr(d)
    idx = d ~= 0;
    g = zeros(size(d));
    g(idx) = d(idx).^(-.5);
    %G = diag(dinv);
end

function y = Afun(x,tflag,A,g1,g2,taus)
    if strcmp(tflag,'notransp')
        y = g1.*(A*(g2.*x)) + taus*g1*(g2'*x);
    else
        y = g2.*(A'*(g1.*x)) + taus*g2*(g1'*x);
    end
end


% function y = Ax(x,A,g1,g2,taus)
%      y = g1.*(A*(g2.*x)) + taus*g1*(g2'*x);
% end
