function F = scfilter(radius, rn, tn, blrate)
%Usage:
%   F = scfilter(rb, tb, blrate)
%Input:
%   radius: radius of the filter
%   rn:     number of radius bins
%   tn:     number of theta (angle) bins
%   blrate: blurring percentage, which indicates how much an angle bin is
%           larger than its nomial size. The nominal size is 2PI / tn.
%Output:
%   F, the shape context filter, a 2*radius+1 by rb*tb 3D matrix. Each page
%      is basicly fan sector like one's embedded on rectangle of zeros.

rx2n1 = radius * 2 + 1; % some constants. Store them to avoid repeated calc

F = zeros( rx2n1, rx2n1, rn * tn); % allocate memory

logr = log( radius + 1);
dlogr = logr / rn; % delta of log radius
logrs = 0 : dlogr : logr; 
r2b = exp( 2 * logrs( 1 : end - 1));
r2e = exp( 2 * logrs( 2 : end ));
% logr2b and e is the beginning and ending log square radius of each radius
% bin

da = 2 * pi / tn; % da is the increment of each theta bin
a = -pi : da : pi;% from -pi to pi
ab = a(1 : end - 1 )- da * blrate*0.5; 
    % ab now is the beginning angle of each blurred theta bin.           
ae = ab + da* (1 + blrate); 
    % ae is the ending angle of each blurred thetat bin.
   
[X Y] = meshgrid( -radius: radius, -radius: radius);
R2 = X.^2 + Y.^2;
theta = atan2(Y, X);

for j = 1 : tn
    for i = 1 : rn
        f = (R2 >= r2b(i) ) .* (R2 < r2e(i)) .* ...
            ( theta >= ab(j) ) .* (theta < ae(j) ); 
        F(:,:, i + (j - 1)*rn ) = f; % column major
    end
end