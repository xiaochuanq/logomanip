function R = corner_detect(I, Sigma)
% corner feature detector using Harris Algroithm
% I is the input in the format of a double grayscale image. The intensity
% of each pixel should range from 0.0 to 255.0
% R is a double matrix contains corner response for each pixel. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A lazy solution:
% R = cornermetric(I, 'Harris');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our implementaiton, almost the same to matlab's method. but we use
% different filter settings and response measures.

sigma = 1.5;
if(nargin >=2 )
    sigma = Sigma;
end
% compute gradients
Ix = imfilter(I,[-1 0 1] ,'replicate','same','conv');
Iy = imfilter(I,[-1 0 1]','replicate','same','conv');

% only use valid gradients
Ix = Ix(2:end-1,2:end-1);
Iy = Iy(2:end-1,2:end-1);

% window function
hsize = max(5, 2* floor(1.5*sigma)+1);
w = fspecial('gaussian', hsize, sigma);

IxIx = imfilter(Ix .* Ix, w, 'replicate', 'full', 'conv');
IxIy = imfilter(Ix .* Iy, w, 'replicate', 'full', 'conv');
IyIy = imfilter(Iy .* Iy, w, 'replicate', 'full', 'conv');

% clip to image size
cut = (hsize-1)/2;
IxIx = IxIx( cut : end - cut + 1, cut : end - cut + 1);
IxIy = IxIy( cut : end - cut + 1, cut : end - cut + 1);
IyIy = IyIy( cut : end - cut + 1, cut : end - cut + 1);

D =  IxIx .* IyIy - 2 * IxIy .* IxIy; % determinant
T =  IxIx + IyIy + eps; % trace

R = D ./ T;