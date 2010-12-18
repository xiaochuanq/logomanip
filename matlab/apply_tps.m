
function [X Y] = apply_tps(cpx, cpy, Wx, Wy, x, y)
% Input
%   Wx, Wy : TPS coefficients. They are m-by-1 
%   cpx, cpy, cooridnates of control points. They are mc-by-1
%   x : the column coords vector, m*1, in the source image
%   y : the column coords vector, m*1, in the source image
% Output
%   X : the column coords vector, m*1, in the destination image
%   Y : the column coords vector, m*1, in the destination image
R = pdist2([x,y], [cpx, cpy], 'euclidean');
% R is m-by-mc if [x y] is m-by-2 and [cpx, cpy] is m-by-2
% each row of D contains [xi,yi]'s distance to each control points
U = R.*R.*log(R + eps); % in most literature U(r) = r*r*ln(r)
K = [ ones(length(x), 1), x, y,  U];
% Attn: W's are in the form of [a1 ax, ay, w]'
X = K * Wx;
Y = K * Wy;
