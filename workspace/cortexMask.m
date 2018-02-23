function mask = cortexMask( level, angle, I )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

imageSizeX = size(I,2);
imageSizeY = size(I,1);
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

centerX = imageSizeX/2;
centerY = imageSizeY/2;

if imageSizeX > imageSizeY
    radius = imageSizeY;
else
    radius = imageSizeX;
end
outradius = radius;
inradius = radius;
switch level
    case 1
        outradius = radius/32;
        inradius = 0;
    case 2        
        outradius = radius/16;
        inradius = radius/32;
    case 3
        outradius = radius/8;
        inradius = radius/16;
    case 4
        outradius = radius/4;
        inradius = radius/8;
    case 5
        outradius = radius/2;
        inradius = radius/4;
    case 6
        outradius = radius/2;
        inradius = radius/2;
end

outCicle = (rowsInImage - centerY).^2 ./ outradius^2 ...
    + (columnsInImage - centerX).^2 ./ outradius^2 >=1;

inCicle = (rowsInImage - centerY).^2 ./ inradius^2 ...
    + (columnsInImage - centerX).^2 ./ inradius^2 >=1;

if level == 1
    mask = 1 - outCicle;
elseif level == 6
    mask = 1 - inCicle; 
else    
    mask = outCicle + (1-inCicle);
    mask = 1 - mask;
end

rectX=[centerX+outradius, centerX+outradius, centerX, centerX-outradius, centerX-outradius];
rectY=[centerY-outradius/3, centerY+outradius/3, centerY, centerY-outradius/3, centerY+outradius/3];
switch angle
    case '0'
        rectX=[centerX+outradius, centerX+outradius, centerX, centerX-outradius, centerX-outradius];
        rectY=[centerY-outradius/3, centerY+outradius/3, centerY, centerY-outradius/3, centerY+outradius/3];
    case '30'
        rectX=[centerX+outradius, centerX+outradius, centerX, centerX-outradius, centerX-outradius];
        rectY=[centerY-outradius/3, centerY-outradius, centerY, centerY+outradius, centerY+outradius/3];
    case '60'
        rectX=[centerX+outradius/3, centerX+outradius, centerX, centerX-outradius, centerX-outradius/3];
        rectY=[centerY-outradius, centerY-outradius, centerY, centerY+outradius, centerY+outradius];
    case '90'
        rectX=[centerX+outradius/3, centerX-outradius/3, centerX, centerX+outradius/3, centerX-outradius/3];
        rectY=[centerY-outradius, centerY-outradius, centerY, centerY+outradius, centerY+outradius];
    case '120'
        rectX=[centerX-outradius, centerX-outradius/3, centerX, centerX+outradius/3, centerX+outradius];
        rectY=[centerY-outradius, centerY-outradius, centerY, centerY+outradius, centerY+outradius];
    case '150'
        rectX=[centerX-outradius, centerX-outradius, centerX, centerX+outradius, centerX+outradius];
        rectY=[centerY-outradius, centerY-outradius/3, centerY, centerY+outradius/3, centerY+outradius];
end

if level > 1 && level < 6    
    maskangle = 1-roipoly(mask,rectX,rectY);
    mask = regionfill(mask,maskangle)==1;
elseif level == 6
    maskangle = 1-roipoly(mask,rectX,rectY);
    imshow(maskangle);
    mask = 1-regionfill(mask,maskangle)==1;
end
end

