function [y1, x1, y2, x2] = extract_match(m, Y1,X1, Y2, X2)
% m is a n1x1 vector of integers where m(i) points to the index of the
% descriptor in p2 that matches with the descriptor p1(:,i). If no match is
% found for feature i, you should put m(i)=-1. 
% p1 & p2 are obtained at Y1,X1 and Y2,X2.
idx = find(m > 0);
x1 = X1(idx);
y1 = Y1(idx);
x2 = X2(m(idx));
y2 = Y2(m(idx));

