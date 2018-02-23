%% HDR-VDP-2: A calibrated visual metric for visibility and quality predictions in all luminance conditions
% matlab version
%
% -----------------------------------------------------------------
% Documentation
%
% HDRVDP - run HDR-VDP-2 on a pair of images
% HDRVDP_VISUALIZE - produce visualization of the probability maps
% HDRVDP_PIX_PER_DEG - compute an angular resolution required for HDR-VDP-2
%
% -----------------------------------------------------------------
% Example: 
%  
clear;
close all;

%%new
ref = imread('/Users/wejaq/Documents/MATLAB/new/1.jpg');
img = initializeImage('',ref);

params = defaultSaliencyParams(img.size);

[salMap,salData] = makeSaliencyMap(img,params);
wta = initializeWTA(salMap,params);
[wta,thisWinner] = evolveWTA(wta);
wtaMap = emptyMap(img.size(1:2),'Winner Take All');
wtaMap.data = imresize(wta.sm.V, img.size(1:2),'bilinear');

threshold = 0.6;
bw = wtaMap.data/max(abs(wtaMap.data(:)));


input = bw;

figure,
imshow(ref,[]);
[B,L] = bwboundaries(input);
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end


minval = 0.05;
maxval = 0.5;

s=10;
A = zeros(size(wtaMap.data,1), size(wtaMap.data,2), 'double');
ys=size(wtaMap.data,1)/s;
xs=size(wtaMap.data,2)/s;

B = wtaMap.data/max(abs(wtaMap.data(:)));  
B = B>threshold;
m = mean(mean(B));

P=[];
for i=0:s-1
    for j=0:s-1
        thisBoundingbox = floor([i*xs+1,j*ys+1,xs,ys]);
        val = mean(mean(imcrop(B, thisBoundingbox)));
        
        P=[val,P];
        if val > minval && val <= maxval 
            A(thisBoundingbox(2):thisBoundingbox(2)+thisBoundingbox(4),...
                thisBoundingbox(1):thisBoundingbox(1)+thisBoundingbox(3)) = 1;
        end
    end
end

minv = min(P);
maxv = max(P);

figure,
imshow(A,[]);

% img = imread( 'sample/nn.jpg' );

%  % Load reference and test images
% img1 = imread( '../sample/11.bmp' );
% img2 = imread( '../sample/33.bmp' );
%   
% %  img1 = img1(71:120, 51:100, :);
% %  img2 = img2(71:120, 51:100, :);
%   
% T = double(img1)/2^8; % images must be normalized 0-1
% R = double(img2)/2^8;
% 
% T1 = rgb2gray(img1);
% R1 = rgb2gray(img2);
% ppd = hdrvdp_pix_per_deg( 21, [size(0,1) size(0,2)], 1 );
% 
% [id,map] = cwssim_index(T1,R1,1,4,0,0);
% figure,imshow(map,[]);
% cwssim = cwssim_index(T1, R1,6,16,0,0);
% res1 = VDM_IDX( img1,img2, 30);
% res1 = hdrvdp_ssim(img1, img2, 'sRGB-display', ppd, [], s, min, max, 'cw'
% map = img1;
% s = size(res.dmap,2);
% % map = cat(3, ssimmap, ssimmap, ssimmap);
% map = double(map)/2^8;
% for i=1:s
%     thisBoundingbox = floor(res.patchPoints(:,i));
%     map(thisBoundingbox(2):(thisBoundingbox(2) + thisBoundingbox(4)),...
%        thisBoundingbox(1):(thisBoundingbox(1) + thisBoundingbox(3)),:) = res.dmap{i};
% end

% figure,imshow(hdrvdp_visualize(res1.dmap));

% [ssimval, ssimmap] = ssim(T1,R1);

% A = zeros(size(img1,1), size(img1,2), 'double');
% 
% % figure,
% % imshow(img1, []);
% tic
% s = 10;
% 
% ys=size(img1,1)/s;
% xs=size(img1,2)/s;
% res=[];
% for i=0:s-1
%     for j=0:s-1
%         thisBlobsBoundingBox = [i*xs,j*ys,xs,ys];
%         if thisBlobsBoundingBox(1)+thisBlobsBoundingBox(3) > size(img1,2)
%             thisBlobsBoundingBox(3)=thisBlobsBoundingBox(3)-5;
%         end
% 
%         if thisBlobsBoundingBox(2)+thisBlobsBoundingBox(4) > size(img1,1)
%             thisBlobsBoundingBox(4)=thisBlobsBoundingBox(4)-5;
%         end
% 
%         if thisBlobsBoundingBox(1) == 0
%             thisBlobsBoundingBox(1)=thisBlobsBoundingBox(1)+3;
%         end
% 
%         if thisBlobsBoundingBox(2) == 0
%             thisBlobsBoundingBox(2)=thisBlobsBoundingBox(2)+3;
%         end
%         
%         patch1 = imcrop(img1, thisBlobsBoundingBox);
%         patch2 = imcrop(img2, thisBlobsBoundingBox);
%         T1 = rgb2gray(patch1);
%         R1 = rgb2gray(patch2);
%         ssimval = ssim(T1,R1);
% %         ssimval = cwssim_index(T1,R1,2,4,0,0);
% %         res(i+1,j+1)=1-ssimval;
%         if 1-ssimval>0.1          
%             A(thisBlobsBoundingBox(2):thisBlobsBoundingBox(2)+thisBlobsBoundingBox(4),...
%        thisBlobsBoundingBox(1):thisBlobsBoundingBox(1)+thisBlobsBoundingBox(3)) = 1;
%         end
%     end
% end
% % 
% % % % CC = imread('1.png');
% % % % c = rgb2gray(CC);
% % % % bw = c < 100;
% % % % bw = 1-bw;
% % % 
% labeledImage = bwlabel(A, 8); 
% % % % A= A.*255;
% Measurements = regionprops(labeledImage, 'Centroid','BoundingBox');
% % toc
% % centroids = cat(1, Measurements.Centroid);
% % hold on
% % plot(centroids(:,1),centroids(:,2), 'b*')
% % hold off
% numberOfBlobs = size(Measurements, 1);
% for k = 1 : numberOfBlobs
%     thisBoundingbox = floor(Measurements(k).BoundingBox);
% 
%     A(thisBoundingbox(2):thisBoundingbox(2)+thisBoundingbox(4),...
%        thisBoundingbox(1):thisBoundingbox(1)+thisBoundingbox(3)) = 1;
% end
% toc
% figure,
% imshow(A, []);

% binaryImage = T1 > 80;
% binaryImage = 1 - binaryImage;
% labeledImage = bwlabel(binaryImage, 8); 
% 
% % mask = binaryImage;
% % 
% % W = graydiffweight(T1, mask, 'GrayDifferenceCutoff', 25);
% % thresh = 0.01;
% % BW = imsegfmm(W, mask, thresh);
% % figure, imshow(BW)
% % title('Segmented Image')
% % cc = bwconncomp(labeledImage, 4)
% % labeled = labelmatrix(cc);
% % RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
% % imshow(RGB_label)
% % binaryImage = 1 - binaryImage;
% % binaryImage = imfill(binaryImage, 'holes');
% 
% figure,
% imshow(labeledImage, []);
% 
% Measurements = regionprops(labeledImage, T1, 'Area', 'BoundingBox');
% numberOfBlobs = size(Measurements, 1);
% for k = 1 : numberOfBlobs
%     if  Measurements(k).Area>400 && Measurements(k).Area<50000
%         thisBlobsBoundingBox = Measurements(k).BoundingBox;
% 
%         hold on
%         rectangle('Position',thisBlobsBoundingBox,'FaceColor','r');
%         hold off
%     end
% end


% ss=1-ssimmap;
% sorted = sort(ss(:),'descend');
% top = sorted(1:1000);
% [~,ia,~] = intersect(ss(:),top(:));
% [r,c]=ind2sub(size(ss),ia);
% 
% Reg = {};
% Point = [];
% meanval=mean2(ss)
% for i=1:length(top)
%     isexist = any(cellfun(@(x) isequal(x, [floor(r(i)/s),floor(c(i)/s)]), Reg));
%     if isexist==0
%         Reg{i}=[floor(r(i)/s),floor(c(i)/s)];
%         Point = [[floor(r(i)/s),floor(c(i)/s)];Point];
%     end
% end
% Point= Point.*s;
% for i = 1:length(r)
%     ss(r(i),c(i)) = 1;
% end
% for i = 1:length(Point)
%     x = Point(i,1)+s;
%     if x>size(ss,1)
%         continue;
%     end
%     y=Point(i,2)+s;
%     if y>size(ss,2)
%         continue;
%     end
%     if mean2(ss(Point(i,1)+1:x, Point(i,2)+1:y))<meanval
%         continue;
%     end
%     ss(Point(i,1)+1:x, Point(i,2)+1:y) = 1;  
% end
% 
% imshow(ss);
% figure, imshow(ssimmap,[]);
% title(sprintf('ssim Index Map - Mean ssim Value is %0.4f',ssimval));

% %
 % Compute pixels per degree for the viewing conditions

 % Run hdrvdp
% tic;
% 
% res = hdrvdp( T, R, 'sRGB-display', ppd );
% b = toc;
% figure,
% imshow( hdrvdp_visualize( res.P_map) );
% 
% s = 10;
% tic;
% res1 = hdrvdp_saliencymap( img1, img2, 'sRGB-display', 20, [], s, 0.05, 0.5, 'normal');
% c = toc;
% if size(res1,1)>0
%     figure,imshow( hdrvdp_visualize_ssim(img1, res1.P_map, res1.patchPoints));
% end
% 
% gray = rgb2gray(img);
% figure,
% imshow(gray);
% F = fftshift(fft2(gray));
% F2 = log(abs(F));
% figure,
% imshow(F2,[]);
% 
% % Create a logical image of an ellipse with specified
% % semi-major and semi-minor axes, center, and image size.
% % First create the image.
% % imageSizeX = size(F,1);
% % imageSizeY = size(F,2);
% % [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% % % Next create the ellipse in the image.
% % centerX = imageSizeX/2;
% % centerY = imageSizeY/2;
% % radius = imageSizeX/32;
% % 
% % ellipsePixels = (rowsInImage - centerY).^2 ./ radius^2 ...
% %     + (columnsInImage - centerX).^2 ./ radius^2 >=1;
% % 
% % ellipsePixels = 1-ellipsePixels;
% 
% % rectX0=[centerX+radius, centerX+radius, centerX, centerX-radius, centerX-radius];
% % rectY0=[centerY-radius/3, centerY+radius/3, centerY, centerY-radius/3, centerY+radius/3];
% % 
% % rectX30=[centerX+radius, centerX+radius, centerX, centerX-radius, centerX-radius];
% % rectY30=[centerY-radius/3, centerY-radius, centerY, centerY+radius, centerY+radius/3];
% % 
% % rectX60=[centerX+radius/3, centerX+radius, centerX, centerX-radius, centerX-radius/3];
% % rectY60=[centerY-radius, centerY-radius, centerY, centerY+radius, centerY+radius];
% % 
% % rectX90=[centerX+radius/3, centerX-radius/3, centerX, centerX+radius/3, centerX-radius/3];
% % rectY90=[centerY-radius, centerY-radius, centerY, centerY+radius, centerY+radius];
% % 
% % rectX120=[centerX-radius, centerX-radius/3, centerX, centerX+radius/3, centerX+radius];
% % rectY120=[centerY-radius, centerY-radius, centerY, centerY+radius, centerY+radius];
% % 
% % rectX150=[centerX-radius, centerX-radius, centerX, centerX+radius, centerX+radius];
% % rectY150=[centerY-radius, centerY-radius/3, centerY, centerY+radius/3, centerY+radius];
% % 
% % mask = 1-roipoly(ellipsePixels,rectX0,rectY0);
% 
% mask = cortexMask(1, '30', F);
% figure,
% imshow(mask);
% % figure,
% % imshow(mask);
% 
% % ellipsePixels = regionfill(ellipsePixels,mask)==1;
% % figure, imshow(ellipsePixels)
% 
% % F=F.*mask;
% % amplitudeImage2 = log(abs(F));
% % ellipsePixels is a 2D "logical" array.
% % Now, display it.
% % imshow(amplitudeImage2,[]);
% % title('Binary image of a ellipse', 'FontSize', 20);
% 
% filteredImage = ifft2(fftshift(F));
% amplitudeImage3 = abs(filteredImage);
% figure,
% imshow(amplitudeImage3,[]);
