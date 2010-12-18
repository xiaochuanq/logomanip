function [M inlier_ind] = ransac_est_model(x1, y1, x2, y2, threshhold, max_iteration)
% function [M inlier_ind] = ransac_est_model(y1, x1, y2, x2, threshhold, max_iteration)
% Input: 
%  x1, y1, the coordinates of source points
%  x2, y2, the coordinates of destination points
%  threshhold: ransac inlier threshold
%  max_iteration, straightforwad, optional and default is 500 for each
%  model. If homography model is good enough, only 500 iterations will be
%  done.
%
% Output:
%  M, the model matrix. If M is 3x3, it is Homography. if it is Nx2
%  it is a TPS model.
%  inlier_ind, the indices of corresponding inlier set of the estimation.
%
% Note: 
%  Source and destination are defined according to the coordinates' mappling
%  directions. For each model, we will find the mapped position of points in
%  image1 (x1, y1) in image2 ([x2 y2] = f(x1, y1)): map from 1 to 2
%  x2, y2 = H * (x1, y1)    ----- f is H * 
%  x2, y2 = TPS (x1, y1)  ----- f is TPS
%
%  However, the pixel values will be transferred reversely, i.e., from the
%  destination to the source. 
%  Later on we can interpolate the pixel ligth/color value in image2 using
%  [x2, y2] = f(x1, y1), and sample back the value for points x1, y1 in
%  image1

thresh = 6; % an egineering guess
if(nargin == 5)
    thresh = threshhold;
end
max_iter = 500;
if(nargin == 6)
    max_iter = max_iteration;
end
n = length(y1);
inlier_ind = [];

if( length(x2)  < 4) 
    disp('Error. Too few points to estimate model');
    return;
end
% firstly we estimate homography. always x2 = H * x1
for iter = 1 : max_iter
    sample = randsample(n, 4); % Task3 step (a)
    H = est_homography(x2(sample), y2(sample), x1(sample), y1(sample)); % step(b)
    [x2hat y2hat] = apply_homography(H, x1, y1); % predict the destinateions
    % x2predict y2predict = H*x, H*y
    % find inliers using threshhold
    new_inlier_ind = find_inlier(x2hat, y2hat, x2, y2, thresh);
    if( length(inlier_ind) < length(new_inlier_ind) );
        inlier_ind = new_inlier_ind; % step e, keep the largest set of inliers
    end
end
% reestimate homography model using the largest set of inliers.
M = est_homography(x2(inlier_ind), y2(inlier_ind), x1(inlier_ind), y1(inlier_ind));

homo_ratio = length(inlier_ind) / n;
if( homo_ratio >= 0.9) % if homography estimation is good enough
    return;            % we won't bother to do the TPS estimate
end

if( length(x2) < 5)
    return;
end

for iter = 1 : max_iter
    sample = randsample(n, 6); % or sample = randsample(n, 6);
    % estimate tps coefficients for x and y respectively
    % in est_tps, we build kernel using ctrl points in image1 and find
    % their transform coefficients W's which map pts in 1(src) to 2(dest).
    W1to2x = est_tps( [x1(sample), y1(sample)], x2(sample));
    W1to2y = est_tps( [x1(sample), y1(sample)], y2(sample));
    %predict x2, y2
    [x2hat y2hat] = apply_tps(x1(sample), y1(sample), W1to2x, W1to2y, x1, y1);    
    % find inliers using threshhold
    new_inlier_ind = find_inlier(x2hat, y2hat, x2, y2, thresh);    
    if( length(inlier_ind) < length(new_inlier_ind) );
        inlier_ind = new_inlier_ind; % step e, keep the largest set of inliers
    end
end
% reestimate homography model using the largest set of inliers.
Wx = est_tps( [x1(inlier_ind), y1(inlier_ind)], x2(inlier_ind));
Wy = est_tps( [x1(inlier_ind), y1(inlier_ind)], y2(inlier_ind));

tps_ratio = length(inlier_ind) / n;
if( tps_ratio >= min(1.0,  homo_ratio * 1.00001) )
    % choose TPS only if it is significantly better than homography model
    M = [Wx, Wy];
end


%% Auxiliary facility functions
function inlier_ind = find_inlier( xhat, yhat, x, y, thresh )
SSD = (xhat - x).^2 + (yhat - y).^2; % step (c)
inlier_ind = find( SSD <= thresh );





