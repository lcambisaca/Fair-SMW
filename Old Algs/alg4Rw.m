
function [clusterLabels,t] = alg4Rw(W, D, F, k) %RW-SWM
%INPUT:
%   W ... (weighted) adjacency matrix of size n x n
%   D ... degree matrix of W
%   F ... group membership matrix G of size n x (h-1)
%   k ... number of clusters
%
%OUTPUT:
% clusterLabels ... vector of length n comprising the cluster label for each
%                  data point
% t ... CPU time of eigs
%-----------------------------------------------------------------------------%

n = size(W, 1);
G = D\W; % Inverse of D always works
M = G*F;
U = G  - M * ((F'*M)\(F'*G'));

tic
[X, vals] = eigs(U, k, 'lr','MaxIterations',1000,'SubspaceDimension', 4*k);
t = toc;
clusterLabels = kmeans(X,k,'Replicates',10, 'MaxIter',500);
end
