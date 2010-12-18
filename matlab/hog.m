function Feature = hog(Img,Pty,Ptx)
%Im : double image/or double matrixs.
%     
%
%Feature: m*n x Feature_siz matrix , Each row is the feature representation
%           of the block 
%
% Recursively adds anything under the  directory to the path.

%fprintf('Adding all subdirectories of current directory to path.\n');
%addpath(genpath(pwd));

%%
%Sizes 
[Row_s Col_s] = size(Img(:,:,1));
Nbins = 9;
Cell_siz  = 6;
Block_siz = 3;
%Block size in pixels
Bl_Pix_siz = Cell_siz*Block_siz;
%Feature size in block 
Feature_siz = Nbins*Block_siz*Block_siz;

%%
%%RGB space
% Gradients kernels
x = [-1 1];
y = x';
%Gradients
Ix = imfilter(Img,x,'same','conv');
Iy = imfilter(Img,y,'same','conv');
Ix(:,end,:) = 0;
Iy(end,:,:) = 0;
%Magnitude and Angles 
Im_Mag = (Ix.*Ix + Iy.*Iy).^0.5;
Im_Angle = atan2(Iy,Ix);
%Signed angles to unsigned angles
%Im_Angle(Im_Angle < 0) = Im_Angle(Im_Angle < 0) + pi;
 %imagesc(uint8(Im_Mag));
%%
%Dominant Orientation
V_x = cos(Im_Angle);
V_y = sin(Im_Angle);
V_x = V_x.*Im_Mag;
V_y = V_y.*Im_Mag;

%%
%Image HoG secroptor basedo on Feature points with dominant orientation
% offset_x = [-13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13];
% offset_x = kron(offset_x,ones(27,1));
% offset_y = offset_x';
% 
% FeaPt_No = numel(Pty);
% 
% patch_x = zeros(27,27,FeaPt_No);
% patch_y = zeros(27,27,FeaPt_No);
% 
% Orient_x = zeros(27*27,1);
% Orient_y = zeros(27*27,1);
% 
% patchA = zeros(27,27);
% patchM = zeros(27,27);
% 
% Patch_Ang = zeros(Bl_Pix_siz*Bl_Pix_siz,FeaPt_No);
% Patch_Mag = zeros(Bl_Pix_siz*Bl_Pix_siz,FeaPt_No);
% 
% 
% for i = 1:1:FeaPt_No
%     patch_x(:,:,i) = Ptx(i)+ offset_x;
%     patch_y(:,:,i) = Pty(i)+ offset_y;  
% end
% 
% patch_x(patch_x < 1) = 1;
% patch_y(patch_y < 1) = 1;
% patch_x(patch_x > size(Img,2)) = size(Img,2);
% patch_y(patch_y > size(Img,1)) = size(Img,1);
% 
% theta = 0;
% for i = 1:1:FeaPt_No
%     ind = sub2ind(size(Img),patch_y(:,:,i),patch_x(:,:,i));
%     ind = reshape(ind,27*27,1);    
%     Orient_x(:) = V_x(ind);
%     Orient_y(:) = V_y(ind);
%     theta = atan2(sum(Orient_y),sum(Orient_x));
%     
%     patchA = reshape(Im_Angle(ind),27,27);
%     patchM = reshape(Im_Mag(ind),27,27);
%     
%     patchA = imrotate(patchA,-theta/pi*180,'crop')-theta;
%     patchM = imrotate(patchM,-theta/pi*180,'crop');
%     
%     patchA(patchA< -pi) = patchA(patchA < -pi) + 2*pi;
%     patchA(patchA> pi) = patchA(patchA > pi) - 2*pi;
%     
%     Patch_Ang(:,i) = reshape(imcrop(patchA,[5,5,17,17]),18*18,1);
%     Patch_Mag(:,i) = reshape(imcrop(patchM,[5,5,17,17]),18*18,1);
% end
% Patch_Ang(Patch_Ang < 0) = Patch_Ang(Patch_Ang < 0) + pi;

%%
%Image HoG descriptor based on Feature points. 
offset_x = [-9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
offset_x = kron(offset_x,ones(18,1));
offset_y = offset_x';

FeaPt_No = numel(Pty);

patch_x = zeros(Bl_Pix_siz,Bl_Pix_siz,FeaPt_No);
patch_y = zeros(Bl_Pix_siz,Bl_Pix_siz,FeaPt_No);
Patch_Ang = zeros(Bl_Pix_siz*Bl_Pix_siz,FeaPt_No);
Patch_Mag = zeros(Bl_Pix_siz*Bl_Pix_siz,FeaPt_No);

for i = 1:1:FeaPt_No
    patch_x(:,:,i) = Ptx(i)+ offset_x;
    patch_y(:,:,i) = Pty(i)+ offset_y;  
end
patch_x(patch_x < 1) = 1;
patch_y(patch_y < 1) = 1;
patch_x(patch_x > size(Img,2)) = size(Img,2);
patch_y(patch_y > size(Img,1)) = size(Img,1);

for i = 1:1:FeaPt_No
    ind = sub2ind(size(Img),patch_y(:,:,i),patch_x(:,:,i));
    ind = reshape(ind,Bl_Pix_siz*Bl_Pix_siz,1);    
    Patch_Mag(:,i) = Im_Mag(ind);
    Patch_Ang(:,i) = Im_Angle(ind);
end
Patch_Ang(Patch_Ang < 0) = Patch_Ang(Patch_Ang < 0) + pi;


%%
%Gray-scale space
% imagesc(uint8(Im));
% Im = rgb2gray(Im);
% x = [-1 1];
% y = x';
% Ix = imfilter(Im,x,'same','conv');
% Iy = imfilter(Im,y,'same','conv');
% Im_Mag = (Ix.*Ix + Iy.*Iy).^0.5;
% Im_Angle = atan2(Iy,Ix);
% 
% Im_Angle(Im_Angle < 0) = Im_Angle(Im_Angle < 0) + pi;
%imagesc(uint8(Ixy));
%quiver(Ix,Iy);

%%
%Feature matrix 
Feature = zeros(FeaPt_No,Feature_siz);

%%



%Block and Cell allocation
Block_Ang = zeros(Cell_siz*Block_siz,Cell_siz*Block_siz);
Block_Mag = zeros(Cell_siz*Block_siz,Cell_siz*Block_siz);
Cell_Ang = zeros(Cell_siz,Cell_siz);
Cell_Mag = zeros(Cell_siz,Cell_siz);


%%
%Don't want it 
% for i = 0:1:m-1
%     for j = 0:1:n-1
%         %Block Magnitude and Block Angles
%         Block_Ang = Im_Angle(i*Bl_Pix_siz+1:(i+1)*Bl_Pix_siz,j*Bl_Pix_siz+1:(j+1)*Bl_Pix_siz,:);
%         Block_Mag = Im_Mag(i*Bl_Pix_siz+1:(i+1)*Bl_Pix_siz,j*Bl_Pix_siz+1:(j+1)*Bl_Pix_siz,:);
%         %%
%         for p = 0:1:Block_siz - 1
%             for q = 0:1:Block_siz - 1
%                 %Cell Magnitude and Cell Angles
%                 Cell_Ang = Block_Ang(p*Cell_siz+1:(p+1)*Cell_siz,q*Cell_siz+1:(q+1)*Cell_siz,:);
%                 Cell_Mag = Block_Mag(p*Cell_siz+1:(p+1)*Cell_siz,q*Cell_siz+1:(q+1)*Cell_siz,:);
%                 %%
%                 %Weighted votes for bins
%                 for w = 0:1:8
%                     Feature(i*m+j+1,(p*Block_siz+q)*Nbins+1+w) = sum(sum(sum(Cell_Mag.*(Cell_Ang>= (w*pi/9)).*(Cell_Ang < ((w+1)*pi/9)))));              
%                 end
%             end
%         end
%         %Normalization 
%         Feature(i*m+j+1,:) = Feature(i*m+j+1,:)/(norm(Feature(i*m+j+1,:),2) + 0.02);
%     end
% end


%%
for i = 1:1:FeaPt_No
    
    %Block Magnitude and Block Angles
    Block_Ang = reshape(Patch_Ang(:,i),Bl_Pix_siz,Bl_Pix_siz);
    Block_Mag = reshape(Patch_Mag(:,i),Bl_Pix_siz,Bl_Pix_siz);
    %%
    for p = 0:1:Block_siz - 1
        for q = 0:1:Block_siz - 1
            %Cell Magnitude and Cell Angles
            Cell_Ang = Block_Ang(p*Cell_siz+1:(p+1)*Cell_siz,q*Cell_siz+1:(q+1)*Cell_siz,:);
            Cell_Mag = Block_Mag(p*Cell_siz+1:(p+1)*Cell_siz,q*Cell_siz+1:(q+1)*Cell_siz,:);
            %%
            %Weighted votes for bins
            for w = 0:1:8
                Feature(i,(p*Block_siz+q)*Nbins+1+w) = sum(sum(sum(Cell_Mag.*(Cell_Ang>= (w*pi/9)).*(Cell_Ang < ((w+1)*pi/9)))));
            end
        end
    end
    %Normalization
    Feature(i,:) = Feature(i,:)/(norm(Feature(i,:),2) + 0.02);
end

 end
