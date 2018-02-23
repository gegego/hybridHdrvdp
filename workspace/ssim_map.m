function [ map,maxv,minv,avgv ] = ssim_map( ref, test, s, ssim_method )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

map = zeros(size(ref,1)+1, size(ref,2)+1, 'double');
ys=size(ref,1)/s;
xs=size(ref,2)/s;
P=[];
for i=0:s-1
    for j=0:s-1
        thisBoundingbox = floor([i*xs+1,j*ys+1,xs,ys]);
        patch_ref = imcrop(ref, thisBoundingbox);
        patch_test = imcrop(test, thisBoundingbox);
        switch lower(ssim_method)
            case 'normal'
                ssimval = ssim(patch_test,patch_ref);
            case 'cw'
                ssimval = cwssim_index(patch_test,patch_ref,1,4,0,0);
        end      
        P=[1-ssimval,P];
        map(thisBoundingbox(2):thisBoundingbox(2)+thisBoundingbox(4),...
            thisBoundingbox(1):thisBoundingbox(1)+thisBoundingbox(3)) = ssimval;
    end
end
P(isnan(P)) = 0;
maxv = max(P);
minv = min(P);
avgv = 1-mean(P);
end

