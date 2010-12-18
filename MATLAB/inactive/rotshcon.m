function SC = rotshcon(SC, Ix, Iy, nr, nt)
% Usage: SC = rotshacon(SC, Ix, Iy, nr, nt), to get the rotation invariant
% robust shape context feature descriptors
%
% Input:
%  SC: general shape context descriptors. Each row contains a feature desc
%  Ix Iy: gradient of sample points
%  nr, nt: the number of radius bins and theta bins used to generate SC
%
% Output:
%  RSC: rotation invariant shape context

da = 2 * pi / nt;  % delta angle
theta = atan2(Iy, Ix) + pi; % theta from 0 to 2pi
idx = floor( theta./da); % get the index of the major orientation 

m = size( SC, 1 );
for i = 1 : m
    sc = reshape( SC(i, :), nr, nt);
    sc = circshift(sc, [0 -idx(i)]);
    SC(i, :) = sc(:);    
end
