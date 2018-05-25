clear all
close all
%read in image
I = fitsread('../Plots/Experiments/Setup 3/fits files/deg90.fits');
%normalise image
I = I - min(I(:));
% gaussian filter image
I2 = imgaussfilt(I,3);
figure;imagesc(I2);
title('sigma = 3')
%enhance contrast
I3 = imadjust(I2);
figure(); imagesc(I3);
title('sigma = 3')

% binarize image 
bw = imbinarize(I3);
figure;
imagesc(bw);

% get rid of backround binary noise
bw1 = bwareaopen(bw,100);
figure;imagesc(bw1);

% Set the spot properties table up and remove central region.
cc = bwconncomp(bw1,4)
spotsdata = regionprops('table',cc,'Area','Centroid','Eccentricity','Orientation');

spotareas = [spotsdata.Area];
[max_area,idx] = max(spotareas);
centre1 = [spotsdata.Centroid(idx,2) spotsdata.Centroid(idx,1)];
spotsdata(idx,:) =[];

% removing the small spots which are likely to throw the centering.
spotareas = [spotsdata.Area];
med = median(spotareas);
standev = std(spotareas);
idxs = find(spotareas < med-2*standev);
spotsdata(idxs,:) =[];


spots = false(size(bw1));
for i = 1:cc.NumObjects
    spots(cc.PixelIdxList{i}) = true;
end
spots(cc.PixelIdxList{idx}) = false;

figure;
imagesc(spots);
hold on
x = 1:640;
Co_effs =  zeros(size(spotsdata,1),2);
for i = 1:size(spotsdata,1)
    theta = spotsdata.Orientation(i);
    x1 = spotsdata.Centroid(i,1);
    y1 = spotsdata.Centroid(i,2);
    m = -tand(theta);
    c = y1+tand(theta)*x1;
    Co_effs(i,:) = [m,c];
    plot(x,-x*tand(theta) + c,'-r')
    hold on
end

numSpots = size(spotsdata,1);
xs = zeros(1,3);
ys = zeros(1,3);
c = combnk(1:numSpots,2);
for i= 1:size(c,1)
    ccell = num2cell(c(i,:));
    [id1,id2] = deal(ccell{:});
    xs(i) = (Co_effs(id2,2)- Co_effs(id1,2))/(Co_effs(id1,1)- Co_effs(id2,1));
    ys(i) = Co_effs(id2,1)*xs(i) + Co_effs(id2,2);
end
centre = Centroid2DGaussian(I);
dispersion_cent = [nanmedian(ys),nanmedian(xs)];
plot(dispersion_cent(2),dispersion_cent(1),'kx');
hold off
