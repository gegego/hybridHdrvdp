function map = vdm_visualize_ssim( background, vdmMap, patchPoints )

%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
map = background;
s = size(vdmMap,2);
% map = cat(3, ssimmap, ssimmap, ssimmap);
map = double(map)/2^8;
map = cat(3, map, map, map);

for i=1:s
    thisBoundingbox = floor(patchPoints(:,i));
    map(thisBoundingbox(2):(thisBoundingbox(2) + thisBoundingbox(4)),...
       thisBoundingbox(1):(thisBoundingbox(1) + thisBoundingbox(3)),:) = hdrvdp_visualize(vdmMap{i});
end

end

