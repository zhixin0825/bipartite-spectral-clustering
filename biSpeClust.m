function [ yh, zh, U, S, V ] = biSpeClust( A, K, L, varargin)
%SPECTRALCLUSTERING Summary of this function goes here
%   Detailed explanation goes here


%%% Optional arguments 
parser = inputParser;
addOptional(parser,'type','US')  % could be 'U' or 'US', defualt is U
addOptional(parser,'kmeans_rep',20)  
addOptional(parser,'rmax', min(K,L))  % maximum number of singular values kept
addOptional(parser,'NORMALIZE', false)  % whether to noramlize rows of the singular vector matrix
addOptional(parser,'RegConst', 1) 
addOptional(parser,'ell1', true) 


parse(parser, varargin{:});
type = parser.Results.type;
kmeans_rep = parser.Results.kmeans_rep;
c = parser.Results.RegConst;
rmax = parser.Results.rmax;
ell1 = parser.Results.ell1;
NORMALIZE = parser.Results.NORMALIZE;

%%
A = regularizeAdj(A,'tau',c,'ell1',ell1);   
[U,S,V] = svds(A,rmax);

if NORMALIZE    
    U = row_l2_normalize(U);
    V = row_l2_normalize(V);
end

% if options.verbose
%     kmopts = statset('Display','iter');
% else
%     kmopts = statset('Display','off'); 
% end
kmeans_options = {'replicates', kmeans_rep, 'onlinephase','off'};%,'start','plus'}; %'start','cluster'};
switch lower(type)
    case 'u' % default
        yh = kmeans(U, K, kmeans_options{:});
        zh = kmeans(V, L, kmeans_options{:});%,'Options',kmopts);
    case 'us'
        yh = kmeans(U*S, K, kmeans_options{:});%,'Options',kmopts);
        zh = kmeans(V*S, L, kmeans_options{:});%,'Options',kmopts);
end

end % biSpectralClustering
% 
% function y = infHandle(x)
% y = x;
% y(isinf(x)) = 0;
% end
% 
% 
% function G = safe_diag_pwr(d)
%     idx = d ~= 0;
%     dinv = zeros(size(d));
%     dinv(idx) = d(idx).^(-.5);
%     G = diag(dinv);
% end
