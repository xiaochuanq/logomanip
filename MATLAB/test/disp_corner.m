function disp_corner(I, Y, X, color)
I = uint8(I);
imshow(I);
hold on;
scatter(X, Y, 'o', color);
hold off;