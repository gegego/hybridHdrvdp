function res = VDM_IDX_Saliency( reference, test, pixels_per_degree, minval, maxval, ssim_method, s, threshold )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ctrspr = 2;    % center spread in min
freqcutcpd = 60/(sqrt(pi)*ctrspr);
            % freq cutoff in cpd
spreadratio = 8;    % surround spread/center spread
ampratio = .685; % surround DC ampl/center DC ampl
consensmax = 114; % contrast sensitivity maximum
maskconthresh = .05; % masked contrast threshold
beta = 2;    % slope

A = zeros(size(reference,1), size(reference,2), 'double');
% load settings
in = load('parameters.mat');
params = in.params;

% make Saliencymap
img_ref = initializeImage('',reference);

if (params.foaSize < 0)
  p = defaultSaliencyParams(img_ref.size);
  params.foaSize = p.foaSize;
end

[salMap,salData] = makeSaliencyMap(img_ref,params);
wta = initializeWTA(salMap,params);
[wta,thisWinner] = evolveWTA(wta);
wtaMap = emptyMap(img_ref.size(1:2),'Winner Take All');
wtaMap.data = imresize(wta.sm.V, img_ref.size(1:2),'bilinear');

%-----------------------------------------------

ys=size(reference,1)/s;
xs=size(reference,2)/s;

norm = wtaMap.data/max(abs(wtaMap.data(:)));
if(threshold>0)
    norm = norm > threshold;
end
% show the result of SaliencyMap
% figure, imshow(norm,[]);

P=[];
for i=0:s-1
    for j=0:s-1
        thisBoundingbox = floor([i*xs+1,j*ys+1,xs,ys]);
        
        val = mean(mean(imcrop(norm, thisBoundingbox)));
            
        P=[val,P];
        if val > minval && val <= maxval           
            A(thisBoundingbox(2):thisBoundingbox(2)+thisBoundingbox(4),...
            thisBoundingbox(1):thisBoundingbox(1)+thisBoundingbox(3)) = 1;
        end
    end
end
% show the result of crop
% figure, imshow(A,[]);

res.maxval = max(P);
res.minval = min(P);

T = (double(reference)-128)/128; % images must be normalized 0-1
R = (double(test)-128)/128;

T1 = rgb2gray(T);
R1 = rgb2gray(R);

rowpixperdeg = pixels_per_degree;
colpixperdeg = pixels_per_degree;

labeledImage = bwlabel(A, 8); 
Measurements = regionprops(labeledImage, 'BoundingBox');
res.same = 1;
for i=1:size(Measurements, 1)
    
    thisBlobsBoundingBox = floor(Measurements(i).BoundingBox);
    
    if thisBlobsBoundingBox(1)+thisBlobsBoundingBox(3) > size(reference,2)
        thisBlobsBoundingBox(3)=thisBlobsBoundingBox(3)-5;
    end
    
    if thisBlobsBoundingBox(2)+thisBlobsBoundingBox(4) > size(reference,1)
        thisBlobsBoundingBox(4)=thisBlobsBoundingBox(4)-5;
    end
    
    if thisBlobsBoundingBox(1) == 0
        thisBlobsBoundingBox(1)=thisBlobsBoundingBox(1)+3;
    end
    
    if thisBlobsBoundingBox(2) == 0
        thisBlobsBoundingBox(2)=thisBlobsBoundingBox(2)+3;
    end
    
    refpatch = imcrop(R1, thisBlobsBoundingBox);
    testpatch = imcrop(T1, thisBlobsBoundingBox);
    
    rows = size(refpatch, 1);
    cols = size(refpatch, 2);

    % Using the image size and resolution parameters, the frequency cut in
    % cycles per degree is converted to row and column frequency cut offs in
    % cycles per image.

    rowfreqcut = rows * freqcutcpd / rowpixperdeg ;
    colfreqcut = cols * freqcutcpd / colpixperdeg ;

    %  The unscaled contrast sensitivity filter  is

    csf = filtdog(rows,cols,rowfreqcut,colfreqcut,ampratio,spreadratio);

    % The unity peak gain contrast sensitivity filter is

    csf = csf/max(max(csf)) ;

    % just noticeable differences between images
    res.patchPoints(:,i)=thisBlobsBoundingBox;
    res.same =0;
    [res.dprime{i},res.dmap{i}] = filtmask(refpatch,testpatch,rowpixperdeg,colpixperdeg,csf,consensmax, ...
        maskconthresh,beta);
end
end

%
% filtmask, the model subroutine
%

function  [dprime,dmap] = filtmask( image1,image2,rowpixperdeg,colpixperdeg,csf, ...
            consensmax,maskconthresh,beta)
        % list of local variables:
% meanlum, con1, con2, fcon1, fcon2, dprime0, dprime

% input images are already contrast images
%meanlum = mean(mean(image1));
%con1 = image1/meanlum - 1;
%con2 = image2/meanlum - 1;

% filtered contrast images
fcon1 = real(ifft2(fft2(image1).*csf));
fcon2 = real(ifft2(fft2(image2).*csf));

% d' without masking
deltaSignal = abs(fcon1 - fcon2);
percentageDifference = deltaSignal ./ fcon2; % Percent by element.

X = grayslice(percentageDifference,16);
% xlswrite('64-128.xlsx',X)
dmap = double(X)/16;
% dmap = ind2rgb(gray2ind(X,255),jet(16));

% figure, imshow(X,jet(16));
% figure, imshow(X,colormap);

dprime0 =consensmax * (sum(sum(abs(fcon1-fcon2).^beta)) / ...
    (rowpixperdeg*colpixperdeg))^(1/beta);

%  background RMS contrast

maskcontrast = sqrt(mean(mean(fcon1.*fcon1)));

%  value returned by the subroutine

dprime = dprime0/sqrt(1 + (maskcontrast/maskconthresh)^2);
end

function result = filtdog(r,c,rf,cf,ar,sr)

result=filtgaus(r,c,rf,cf)-ar*filtgaus(r,c,rf/sr,cf/sr) ;

end
%
% filtgaus: 2D Gaussian filter (outer product of two 1D Gaussians)
%

function  result = filtgaus(r, c, rf, cf)

result=fltgaus1(r,rf)'*fltgaus1(c,cf);
end
%
% fltgaus1: generate a 1D Gaussian low pass filter
%
% n: spatial image length (must be even)
% f: 1/e frequency cutoff in cycles per image

function result = fltgaus1( n, f)
if rem(n,2)==0
    result = [[0:n/2] [1:n/2-1]-n/2];
else
    result = [[0:n/2] [1:n/2]-n/2];
end
result = exp(-(result.*result)/(f*f))  ;
%
end