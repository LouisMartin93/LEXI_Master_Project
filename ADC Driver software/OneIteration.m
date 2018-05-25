function dispersion_cent = OneIteration(videoObject,DarkFrame)
%% Pre set constants
lD2px = 3.2; % pixels per lambda/D
r_search = 5*lD2px; % ROI for finding spot maxima (px)
r_cent = 3*lD2px; % ROI for calculating centroid (px)
centClip = 5; % Lower clip threshold for centroiding (peak flux / centClip)

%% Take first image
image = im2double(take_exposure(videoObject));
image = abs(image- DarkFrame);

addpath('spot extraction noGUI/bin');
% plot image for initial inspection
figure();
imagesc(log10(image)); hold on % for measured data use the following limits [-2.5,max(max(log10(abs(image))))]
colorbar(); hold on

decision(videoObject)
disp('Ready to go')
close;
image = im2double(take_exposure(videoObject));
image = abs(image- DarkFrame);
% plot image after inspection
figure();
imagesc(log10(image)); hold on % for measured data use the following limits [-2.5,max(max(log10(abs(image))))]
colorbar(); hold on

%% Display interactive plot
nModes_click = input('Number of spots to select: ');
figure();
norm_data = (image - min(image(:)))/(max(image(:))-min(image(:)));
set(gca,'YDir','normal');
imagesc(log10(abs(norm_data)),[-2.5,-1]); hold on % scaled display
title(strcat('Nmodes: ',num2str(nModes_click),'- click on all blue copies first')); hold on 
    
% Get user input clicks
[x,y] = ginput(nModes_click*2); % Get user clicks for all modes
yp = y(1:nModes_click);
ym = y(nModes_click+1:end);
xp = x(1:nModes_click);
xm = x(nModes_click+1:end);
spotCoords_userInput = [yp,xp,ym,xm];
close();
%% Calculate the Dispersion Centre

%refine selection
spotCoords_px = refine_coords2(image, spotCoords_userInput, r_search, r_cent, centClip);

%display points
plot(spotCoords_px(:,2),spotCoords_px(:,1),'bx'); hold on % Blue
plot(spotCoords_px(:,4),spotCoords_px(:,3),'rx'); hold on % red

nModes_sense = size(spotCoords_px,1);
%calculate line coeffs
Co_effs =  zeros(nModes_sense,2); % Line coefficients (mode (+/- spots) , m/c)
for i = 1:nModes_click
    xb = spotCoords_px(i,2);
    yb = spotCoords_px(i,1);
    xr = spotCoords_px(i,4);
    yr = spotCoords_px(i,3);
    m = (yb-yr)/(xb-xr);
    c  = yb-m*xb;
    Co_effs(i,:) = [m,c];
    %x = 1:640;
    %y = m*x+c;
    %plot(x,y,'k');
end

xs = zeros(1,3);
ys = zeros(1,3);
c = combnk(1:nModes_click,2);
for i= 1:size(c,1)
ccell = num2cell(c(i,:));
[id1,id2] = deal(ccell{:});
xs(i) = (Co_effs(id2,2)- Co_effs(id1,2))/(Co_effs(id1,1)- Co_effs(id2,1));
ys(i) = Co_effs(id2,1)*xs(i) + Co_effs(id2,2);
end
dispersion_cent = [nanmedian(ys),nanmedian(xs)];
plot(dispersion_cent(2),dispersion_cent(1),'kx');
end
