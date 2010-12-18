function disp_inlier_match(y1, x1, y2, x2, im1rgb, im2rgb, inlier_ind)
X1 = x1(inlier_ind);
Y1 = y1(inlier_ind);
X2 = x2(inlier_ind);
Y2 = y2(inlier_ind);
subplot(1,2,1); disp_corner(im1rgb, Y1, X1); 
subplot(1,2,2); disp_corner(im2rgb, Y2, X2);