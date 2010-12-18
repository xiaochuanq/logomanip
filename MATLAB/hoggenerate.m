%function feature = hoggenerate(Templ_Img)

Img_o =double(rgb2gray(imread('B.jpg')));
Img_o(Img_o == 0) = 1;
Scale = zeros(6,1);
%figure,imagesc(Img_o),colormap(gray);

Gau = fspecial('gaussian', 3, 0.6);

Img = cell(6,9);
%%
for i = 1:1:6;
    Img{i,5} = imresize(Img_o,i*0.1+0.4);
    Scale(i) = i*0.1+0.4;
    %imagesc(Img{i,5}),colormap(gray);
    for j = 1:1:9;
        Img{i,j} = imrotate(Img{i,5},-40+(j-1)*10,'crop');
        Img{i,j}(Img{i,j} == 0) = max(max(Img{i,j}));
        Img{i,j} = imfilter(Img{i,j}, Gau, 'replicate', 'same', 'conv');
    end
    %pause(0.1);
end
%%
R = cell(6,9);
Pty = cell(6,9);
Ptx = cell(6,9);
Feature = cell(6,9);
%
disp('Generating Features for all samples');
t0 = clock();
for i = 1:1:6
    for j = 1:1:9;
        R{i,j} = corner_detect(Img{i,j});
        thre = 0.1*max(max(R{i,j}));
        R{i,j}(R{i,j}<thre) = 0;
        [Pty{i,j} Ptx{i,j} rmax] = anms(R{i,j},500);
        Pty{i,j} = round(Pty{i,j});
        Ptx{i,j} = round(Ptx{i,j});
        Feature{i,j}= HOG(Img{i,j},Pty{i,j},Ptx{i,j});
    end
    eta((i-1)*9+j,54,t0);
end
%%
Centroid = cell(6,9);
for i =1:1:6
    for j = 1:1:9
        Centroid{i,j} = [round(0.5*size(Img{i,j},1)) round(0.5*size(Img{i,j},2))]; 
    end
end
%%
Model.Feature = Feature;
Model.Ptx = Ptx;
Model.Pty = Pty;
Model.Centroid = Centroid;
Model.Img = Img;
Model.Scale = Scale;

save B.mat -struct Model
