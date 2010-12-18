function [m] = feat_match(p1, p2)

%p1 p2 see requirement of knnsearch

[m, dist] = knnsearch( p2, p1, 'K', 2, 'Distance', 'Euclidean', 'NSMethod', 'kdtree');
%Attn: knnsearch(X, Y) find nn's in X for each Y.


ratio = dist(:,1)./dist(:,2); % dist 1 is supposed to be samller than dist 2
max_ratio = max(ratio);
min_ratio = min(ratio);
m = m(:,1);
n = length(m);
h = hist(ratio, 100);
percentage = 0.1; % we want to select the top 10% best matches
threshhold = ( find( cumsum(h) > n * percentage, 1) )/100 * (max_ratio - min_ratio) + min_ratio;
m( ratio > threshhold ) = -1; 

