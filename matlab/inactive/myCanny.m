function [E Ix Iy Imag] = myCanny(I, ratio_nm, ratio)
% I: input image
% poe: rotio of the pixles passed the non maximum suppression can be edge
% rotio: low threshold / high threshold
    [Imag Ix Iy] = myImgrad(I);
    global edge;
    [idxLM cordX cordY dX dY] = non_max_sup(Ix, Iy, Imag, 'linear');% non-max suppression   
    maxMag = max(Imag(idxLM));    
    minMag = min(Imag(idxLM));
    numPassed = length(idxLM);
    counts = hist(Imag(idxLM), 100);
    highthresh = find(cumsum(counts) > (1 - ratio_nm) * numPassed, 1, 'first')/100 * (maxMag - minMag)+minMag;
    lowthresh = ratio * highthresh;
    hysteresis(Imag, idxLM, cordX, cordY, dX, dY, highthresh, lowthresh);
    E = (edge > 2.5);

function [idxLM cordX cordY dX dY] = non_max_sup(Ix, Iy, Imag, interptype)
%input, "derivatives" in different directions: Ix, Iy
%gradients: Imag
%how we interpolate: interptype, a string
%output, the index of local maximums, idxLM
%
%O---(j-1)--------j-----------(j+1)--> X
%|
%i-1   o [case 2) o [case 3) o
%|
%|(case 1]                  [case 4)
%|
%i     o          o          o
%|
%|(case 4]                    [case 1)
%|
%i+1   o (case 3] o (case 2] o
%|
%v
%Y
    [m n] = size(Imag);
    [cordX cordY] = meshgrid(1:n, 1:m);
    dX = ones(m, n);
    dY = ones(m, n);
    cord1X = zeros(m, n);
    cord1Y = zeros(m, n);
    cord2X = zeros(m, n);
    cord2Y = zeros(m, n);
%case 1 and 4:
    idx = find( (Iy >= 0 & Ix > Iy ) | (Iy <= 0 & Ix < Iy) );
    idx = [idx; find( (Iy >0 & Ix <=-Iy) | (Iy < 0 & Ix > -Iy) )];
    dY(idx) = Iy(idx)./Ix(idx);
    cord1X(idx) = cordX(idx) + 1;
    cord1Y(idx) = cordY(idx) + dY(idx);
    cord2X(idx) = cordX(idx) - 1;
    cord2Y(idx) = cordY(idx) - dY(idx);
%case 2 and 3:
    idx = find( (Ix > 0 & Iy >= Ix) | (Ix < 0 & Iy <= Ix) );
    idx = [idx; find( (Ix <= 0 & Iy >-Ix) | (Ix >= 0 & Iy <-Ix) )];
    dX(idx) = Ix(idx)./Iy(idx); 
    cord1X(idx) = cordX(idx) + dX(idx);
    cord1Y(idx) = cordY(idx) + 1;
    cord2X(idx) = cordX(idx) - dX(idx);
    cord2Y(idx) = cordY(idx) - 1;  
% interpolate    
    mag1 = interp2(cordX, cordY, Imag, cord1X, cord1Y, interptype);
    % We are luck here. The interpolations go out of the image boundary
    % did not throw any exceptions by the function interp2
    mag2 = interp2(cordX, cordY, Imag, cord2X, cord2Y, interptype);
 % find local maximum   
    idxLM = find(Imag > mag1 & Imag > mag2);
   

function hysteresis(Imag, idxLM, cordX, cordY, dX, dY, ht, lt)   
    global edge;    
    [m n] = size(Imag);
    edge = zeros(m, n); % 0 for not an edge
    edge( idxLM(Imag(idxLM) >= lt) ) = 1; % 1 for a possible edge
    edge( idxLM(Imag(idxLM) >= ht) ) = 2; % 2 for a strong edge
    for i = 2 : m-1
        for j = 2 : n-1
            if edge(i, j) == 2  % go over all strong edge as starter points
                edge(i, j) = 3; % mark as visitied.                
                dir_trace(i, j, cordX, cordY, dX, dY); % follow this starter point
            end
        end
    end
    
function dir_trace(i, j,cordX, cordY, dX, dY)
    global edge;
    [m n] = size(edge);
    step = [1 -1];%  2 -2 3 -3];
    nextX = cordX(i, j) + dX(i, j) * step;
    nextY = cordY(i, j) + dY(i, j) * step;
    for k = 1 : length(step)
        x = round(nextX(1,k));
        y = round(nextY(1,k));
        if( y > 0 ) && ( y <= m) && ( x > 0 ) && ( x <= n)
            if  edge(y,x) ==1 % || edge(y, x) == 2
                edge(y,x) = 3;
                dir_trace( y, x, cordX, cordY, dX, dY);
            end
        end
    end
    
    
function [grad gx gy] = myImgrad(I,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This Function verifies convolution is indeed assoicative and compute the gradient of the image 
%Author: Xiaochuan Qin, Fan Deng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image I should be a graysacle picture

%Smooth
Gx = [0.05 .25  .4 .25 .05];
Gy = Gx';

%Difference Kernel 
%%
if nargin < 2 
    mode = 'diff';
end
%%
if( strcmp(mode, 'sobel'))
deltax = [1 2 1; 
          0 0 0;
          -1 -2 -1];
else if  strcmp(mode, 'diff')
deltax = [1 -1];
end
deltay = deltax';
gx = conv2( conv2( conv2(I, Gx, 'same'), deltax, 'same'), Gy, 'same');
gy = conv2( conv2( conv2(I, Gy, 'same'), deltay, 'same'), Gx, 'same');
grad = sqrt( gx.*gx + gy.*gy);

end

