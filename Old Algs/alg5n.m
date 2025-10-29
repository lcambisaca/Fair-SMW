
function [clusterLabels , t] = alg5n(W, D, F, k) % W + D Doesnt work
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
G = W + diag(D);
M = G*F;
U = G  - M * ((F'*M)\M');

tic;
% [X, vals] = eigs(U, k, 'lr','MaxIterations',1000,'SubspaceDimension', 4*k);
[X, vals] = eigs(U, k, 'lr' ,'MaxIterations',1000);
% [X, vals] = eigs(@(b) SIAM_Afun(b,G,F,inv1), n , k, 'lr','MaxIterations',1000,'SubspaceDimension', 4*k);
t = toc;
clusterLabels = kmeans(X,k,'Replicates',10, 'MaxIter',500);

end



