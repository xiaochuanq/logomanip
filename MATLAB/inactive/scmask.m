function M = scmask( edge, X, Y, r, thresh)
% generate mask for sampled edge points in [X Y].
% r is the radius of the shape context filter.
% thresh is the 

[m n] = size(edge);
E = zeros( m + r, n + r);
BW = zeros( m + r, n + r);
E(1 + r : 1 + m + r, 1+r, 1+ n +4) = edge(:, :);
D = bwdist( E, 'euclidean');
idx_far = find(D > thresh);
idx_near = find(D <= thresh);
BW(idx_far) = 0;
BW(idx_near) = 1;

num = numel(X);
M = zeros(num * (2*r+1), (2*r+1) );

