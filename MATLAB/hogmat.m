function hogmat(img, scale, rotx, roty, rotz)
% img: input image, in double, one chanel. Before transform, the center of
% the image is placed at (0, 0, 0) conceptually.
%
% scale: a vector of number ranges from (0, inf)
%
% rotx, roty, rotz, the vectors contains rotation angles along x, y, z axis
% respectively. rotx and roty range (0, pi/2), but rotz can range form
% -pi to pi
%
% templates, a cell array. one row for each sacle level, and each column
% for each rotation angle.
% transforms, a cell array contains the tform objects.
%     all tform objects can be used to vote center later on.
% Knowing dx, dy, the vector from [x, y] to [x0, y0]. we can vote the
% center.
% By definition [dx dy] = [x0 y0] - [x y] ( in the original template coord)
% if [x' y'] matches [x y], then a hypothesis is that T*[x y] = [x' y'], (T
% is a tform object). So are T*[x0, y0] = [x0' y0']. 
% T*[dx dy] = T*[x0 y0] - T*[x y]; then [x0' y0'] = T([dx dy] + [x y]) = 
% [x' y'] + T*[dx dy]

[m n] = size(img);
r = (n + 1)/2;
maskimg = zeros(m,n);
center = [ (m+1)/2, (n+1)/2];
[X Y] = meshgrid(1:n, 1:m);
R2 = ( (max(m,n)+1)/2).^2;
maskimg( X.^2 + Y.^2 <= 260*260) = 1;
scale = scale(:);
rotx = rotx(:);
roty = roty(:);
rotz = rotz(:);
fillcolor = max(max(img));

M = length(scale);
nseg = 1 + numel(rotx)*2 + 2*numel(roty) ; % segment length
npos = ( 0 : ( numel(rotz) -1) ) * nseg  + 1; % starting index of each segment
N = nseg * numel( rotz );
%N = n* (numel(crate) + 1);


templates = cell(M, N); % the deformed template logos
transforms = cell(M, N); % and corresponding TForm objects
masks = cell(M, N);
Feature = cell(M, N); % HOG features
Ptx = cell(M, N); % where HOG features are extracted
Pty = cell(M, N);
Dx = cell(M, N); % their relative position to geometric center
Dy = cell(M, N); 
Corner = cell(M, N);
Centroid = cell(M, N);

for i = 1 : M
    scale(i)
    S = eye(3) * scale(i);
    S(3, 3) = 1;  % scaling matrix;
    for rz = 1 : numel( rotz )
        rotz(rz)
        theta = rotz(rz);
        c = cos( theta);
        s = sin( theta);
        Rz = [ c -s 0; s c 0; 0 0 1]; % rotation matrix along z-axis
        RzS = Rz * S;
        tform = maketform('affine', RzS);
        j = npos(rz);
        transforms(i, j ) = {tform};
        templates(i, j) = {imtransform(img, tform, 'bicubic', 'fill', fillcolor)};
   %     masks(i, j) = {imtransform(maskimg, tform, 'bicubic', 'fill', 0)};
        Corner(i, j) = { corner_dect( templates{i, j}  )};
        thresh = 0.1 * max(max(Corner{i, j}));
        Corner{i,j}(Corner{i,j}<thresh) = 0;
        [Pty{i, j} Ptx{i, j} rmax] = anms(Corner{i,j}, 500);
        Ptx{i, j} = round( Ptx{i, j} );
        Pty{i, j} = round( Pty{i, j} );
        Centroid{i,j} =  floor( (size(templates{i, j}) + 1)/2) ;
        Dx{i, j} = Centroid{i, j}(2) - Ptx{i, j};
        Dy{i, j} = Centroid{i, j}(1) - Pty{i, j};
        Feature{i, j} = HOG( templates{i, j}, Pty{i,j}, Ptx{i,j});
        for rx = 1 : numel( rotx)
            theta = rotx(rx);
            Rx = eye(3);
            Rx(2,2) = cos(theta); % rotation matrix along x-axis
            tform = maketform('affine', Rz*Rx  *S);
            j = npos(rz) + rx;
            templates{i,j} = imtransform(img, tform,'bicubic', 'fill', fillcolor);
            transforms( i, j) = {tform};
            masks(i, j ) = {imtransform(maskimg, tform, 'bicubic', 'fill', 0)};
            Corner(i, j) = { corner_dect( templates{i, j}  )};
            thresh = 0.1 * max(max(Corner{i, j}));
            Corner{i,j}(Corner{i,j}<thresh) = 0;
            [Pty{i, j} Ptx{i, j} rmax] = anms(Corner{i, j}, 500);
            Ptx{i, j} = round( Ptx{i, j} );
            Pty{i, j} = round( Pty{i, j} );
            Centroid{i,j} = floor( (size(templates{i, j}) + 1)/2);
            Dx{i, j} = Centroid{i, j}(2) - Ptx{i, j};
            Dy{i, j} = Centroid{i, j}(1) - Pty{i, j};
            Feature{i, j} = HOG( templates{i, j}, Pty{i,j}, Ptx{i,j});           
        end
        for rx = 1 : numel( rotx)
            theta = rotx(rx);
            Rx = eye(3);
            Rx(2,2) = cos(theta); % rotation matrix along x-axis
            tform = maketform('affine', Rx* Rz  *S);
            j = npos(rz) +numel(rotx) + rx;
            templates{i,j} = imtransform(img, tform,'bicubic', 'fill', fillcolor);
            transforms( i, j) = {tform};
            masks(i, j ) = {imtransform(maskimg, tform, 'bicubic', 'fill', 0)};
            Corner(i, j) = { corner_dect( templates{i, j}  )};
            thresh = 0.1 * max(max(Corner{i, j}));
            Corner{i,j}(Corner{i,j}<thresh) = 0;
            [Pty{i, j} Ptx{i, j} rmax] = anms(Corner{i, j}, 500);
            Ptx{i, j} = round( Ptx{i, j} );
            Pty{i, j} = round( Pty{i, j} );
            Centroid{i,j} = floor( (size(templates{i, j}) + 1)/2);
            Dx{i, j} = Centroid{i, j}(2) - Ptx{i, j};
            Dy{i, j} = Centroid{i, j}(1) - Pty{i, j};
            Feature{i, j} = HOG( templates{i, j}, Pty{i,j}, Ptx{i,j});           
        end
        for ry = 1 : numel( roty)
            theta = roty(ry);
            Ry = eye(3);
            Ry(1,1) = cos(theta); % rotation matrix along y-axis
            tform = maketform('affine', Ry* Rz* S);
            j = npos(rz) + numel(rotx)*2 + ry;
            templates(i, j) ={imtransform(img, tform,'bicubic', 'fill', fillcolor)};
            transforms( i, j) = {tform};
            masks(i, j) = {imtransform(maskimg, tform, 'bicubic', 'fill',0)};
            Corner(i, j) = { corner_dect( templates{i, j}  )};
            thresh = 0.1 * max(max(Corner{i, j}));
            Corner{i,j}(Corner{i,j}<thresh) = 0;
            [Pty{i, j} Ptx{i, j} rmax] = anms(Corner{i, j}, 500);
            Ptx{i, j} = round( Ptx{i, j} );
            Pty{i, j} = round( Pty{i, j} );
            Centroid{i,j} = floor( (size(templates{i, j}) + 1)/2);
            Dx{i, j} = Centroid{i, j}(2) - Ptx{i, j};
            Dy{i, j} = Centroid{i, j}(1) - Pty{i, j};
            Feature{i, j} = HOG( templates{i, j}, Pty{i,j}, Ptx{i,j});
        end    
        for ry = 1 : numel( roty)
            theta = roty(ry);
            Ry = eye(3);
            Ry(1,1) = cos(theta); % rotation matrix along y-axis
            tform = maketform('affine', Rz* Ry* S);
            j = npos(rz) + numel(rotx)*2 + numel(roty) + ry;
            templates{i, j} =imtransform(img, tform,'bicubic', 'fill', fillcolor);
            transforms( i, j) = {tform};
            masks(i, j) = {imtransform(maskimg, tform, 'bicubic', 'fill',0)};
            Corner(i, j) = { corner_dect( templates{i, j}  )};
            thresh = 0.1 * max(max(Corner{i, j}));
            Corner{i,j}(Corner{i,j}<thresh) = 0;
            [Pty{i, j} Ptx{i, j} rmax] = anms(Corner{i, j}, 500);
            Ptx{i, j} = round( Ptx{i, j} );
            Pty{i, j} = round( Pty{i, j} );
            Centroid{i,j} = floor( (size(templates{i, j}) + 1)/2);
            Dx{i, j} = Centroid{i, j}(2) - Ptx{i, j};
            Dy{i, j} = Centroid{i, j}(1) - Pty{i, j};
            Feature{i, j} = HOG( templates{i, j}, Pty{i,j}, Ptx{i,j});
        end 
    end 
end

%%
Model.Mask = masks;
Model.Transform = transforms;
Model.Feature = Feature;
Model.Ptx = Ptx;
Model.Pty = Pty;
Model.Dx = Dx;
Model.Dy = Dy;
Model.Centroid = Centroid;
Model.Img = templates;
Model.Scale = scale;

% %% Debug Test
% for i = 1 : M
%     for j = 1:N
%           subplot(M, N, (i-1)*N + j);
% %         size(Model.Img{i,j})
%         imshow(uint8(Model.Img{i,j}));
% %         s = sprintf('logo%d%d.jpg', i, j);
% %         imwrite(uint8(Model.Img{i,j}), s, 'jpg');
%         hold on;
%         plot(Centroid{i, j}(2), Centroid{i,j}(1),'o');
%     end
% end

save B.mat -struct Model

