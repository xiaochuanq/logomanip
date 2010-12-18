function logorematch(imgfile, logofile, clipartfile, modelfile, colfile, outfile, outputfile)

Imorigin = double( imread(imgfile) );
Imlogo = double( rgb2gray( imread( logofile)));
Imclip = double( imread(clipartfile) );

Model = load(modelfile);
column =  load(colfile);
Found = load(outfile);
Found.Col = column.Column;

num = length( Found.Scale );
for i = 1 : num
    s =  Found.Scale(i) ; %  to index Model
    c = Found.Col(i);     %   to index Model
    imtest = Found.Im_out(:,:, i);
    bw = Found.BW(:,:, i);
    imlogo = Model.Img{s, c}(:,:); % no need to scale a logo
    
    % scale first,
    scale = Model.Scale( s );
    imclip = imresize( Imclip, scale);
    imlogo = imresize( Imlogo, scale);
    % and then translate
    theta = -40 + (c - 1)* 10;
    imclip = imrotate( imclip, theta, 'bicubic');
    imlogo = imrotate( imlogo, theta, 'bicubic');
    
    feat_logo = Model.Feature{s, c}; 
    Ylogo = Model.Pty{s, c};
    Xlogo = Model.Ptx{s, c};
    % making the mask. there are many ways to do this
    e = edge(uint8(imlogo),'canny');
    [masky maskx] = ind2sub( size(e), find( e > 0));
    % insert clip art to original image at appropreate position
    Imorigin = addclip(Xlogo, Ylogo, feat_logo, ...
        Imorigin, imclip, imtest, bw, [maskx masky]);
end
imshow(uint8(Imorigin));
imwrite( uint8(Imorigin), outputfile, 'jpg');

