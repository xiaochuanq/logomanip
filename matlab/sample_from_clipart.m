function imgr = sample_from_clipart(imorigin, fpts, bw, imgclip, Model, clipmask)
%input: imorigin: the image to be tested 3 chanel
%       fpts: feature points, in x, y format
%       imgclip: the clip art (3 chanel)
%       clipmask: the mask put on clipart
%output: imgr
%

if( nargin < 5) %default clipmask is the whole area of the clip art
    m = size(imgclip, 1);
    n = size(imgclip, 2);
    clipmask = [1 1; n, 1; n, m; 1, m]; % in (x y) pair, clockwise order
end

% 
% K = convhull(fpts(:,1), fpts(:,2));
% l = min( fpts(K, 1));
% r = max( fpts(K, 1));
% t = min( fpts(K, 2));
% b = max( fpts(K, 2));

idx = find(bw > 0);
[pty ptx] = ind2sub( size(bw), idx);
% K = convhull(ptx, pty);
l = min( ptx(:) );
r = max( ptx(:) );
t = min( pty(:) );
b = max( pty(:) );

[x y] = meshgrid(l:r, t:b);
[m n] = size(x);
imblock1 = imorigin(t:b, l:r, :);

x = x(:);
y = y(:);
[in1 on1]= inpolygon(x, y, [ l r r l l], [ t t b b t ]);
in1 = in1 + on1;

if( size(Model, 2) == 3) % homography
    [X Y] = apply_homography(Model, x, y);
else %TPS
    %fpts = [fpts; l t; r t; r b; l b];
    [X Y] = apply_tps(fpts(:,1), fpts(:,2), Model(:,1), Model(:,2), x, y);
end

k = convhull(clipmask(:,1), clipmask(:,2) ); %convhull(x,y)
[in2 on2 ]= inpolygon( X, Y, clipmask(k, 1), clipmask(k,2));
in2 = in2+on2;

imblock2 = zeros(m, n);

X = reshape(X, m,n);
Y = reshape(Y, m,n);
imblock2(:,:,1) = interp2(imgclip(:,:,1), X, Y, 'cubic', 255);
imblock2(:,:,2) = interp2(imgclip(:,:,2), X, Y, 'cubic', 255);
imblock2(:,:,3) = interp2(imgclip(:,:,3), X, Y, 'cubic', 255);

in1 = reshape(in1, m, n);
in2 = reshape(in2, m, n); 
in = in2;% & in1;

imgr = imorigin;
imgr(t:b, l:r, 1) = (1 - in).* imblock1(:,:,1) + in.* imblock2(:,:,1);
imgr(t:b, l:r, 2) = (1 - in).* imblock1(:,:,2) + in.* imblock2(:,:,2);
imgr(t:b, l:r, 3) = (1 - in).* imblock1(:,:,3) + in.* imblock2(:,:,3);



