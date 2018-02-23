function map = hdrvdp_visualize_ssim( background, hdrvdpMap, patchPoints )

s = size(hdrvdpMap,2);
% map = cat(3, ssimmap, ssimmap, ssimmap);
map = double(background)/2^8;
map = cat(3, map, map, map);
for i=1:s
    thisBoundingbox = floor(patchPoints(:,i));
    if(thisBoundingbox == [0; 0; 0; 0])
        continue;
    end
    map(thisBoundingbox(2):(thisBoundingbox(2) + thisBoundingbox(4)),...
       thisBoundingbox(1):(thisBoundingbox(1) + thisBoundingbox(3)),:) = hdrvdp_visualize( hdrvdpMap{i});
end
end