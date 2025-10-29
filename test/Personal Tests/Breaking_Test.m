clc; clear
format short

n = 4;         % Number of nodes
k = 2;            % Number of clusters
h = 2;            % Number of groups

W = [0 1 0 1; 1 0 1 0 ; 0 1 0 1 ; 1 0 1 0];
% W = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0;];
%W = [0 1 2 3; 1 0 2 3; 2 2 0 3; 3 3 3 0]

% W = [0 1 0 1 0; 1 0 1 1 1; 0 1 0 0 0; 1 1 0 0 1; 0 1 0 1 0]

% W = [0 20 20 20; 20 0 0 0; 20 0 0 0; 20 0 0 0]; % Example adjacency matrix

D = diag(sum(W, 2));

f = [1 0; 0 1 ;0 1; 1 0 ; 0 1];
F = zeros(n, h-1);

for g = 1:(h-1)
    col_sum = sum(f(:, g));
    for h = 1:n
        F(h,g) = f(h,g) - col_sum/n;
    end
end

% clusterLabels2 = alg4(W, D, F, k)
clusterLabels3 = alg3(W, D, F, k)
clusterLabels5 = alg5(W, D, F, k)
clusterLabels6 = alg6Rw(W, D, F, k)
clusterLabels6mod = alg6Sym(W, D, F, k)




