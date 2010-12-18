function D = chi2dist(P, Q)
% Usage: D = chi2dist(P, Q)
% Input: 
%  P, Q, mp-by-descriptor_size and mq-by-descriptor_size feature matrices
%  D, mp-by-mq distance matrix

D = pdist2(P,Q, @mychi2);


function c = mychi2(p, Q)
m = size(Q, 1);
P = repmat(p,[m, 1]);
N = (P - Q).^2;
D =  P + Q + eps;
C = N ./ D;
c = 0.5 * sum( C, 2);
