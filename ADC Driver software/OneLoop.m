%% Script of one loop of ADC
% Assumes that:
% dark image is called MasterDark.
% h1 and h2 are global 
% The ADC angle is known

%% Take first image
image = im2double(take_exposure(vid));
image = abs(image- MasterDark);

%% Run the spot extraction

addpath('spot extraction noGUI/bin');

%input parameters (predefined)
lD2px = 3.2; % pixels per lambda/D
imSize = [640,480]; % pixels (not used)
bandwidth = 145; %nm (not used)
lambda_cent = 560.5; % Central wavelength (nm) (not used)
nModes_click = 8; % Number of modes to manually select
r_search = 5*lD2px; % ROI for finding spot maxima (px)
r_cent = 3*lD2px; % ROI for calculating centroid (px)
centClip = 5; % Lower clip threshold for centroiding (peak flux / centClip)

figure();
%plot image
imagesc(log10(image)); hold on % for measured data use the following limits [-2.5,max(max(log10(abs(image))))]
colorbar(); hold on


figure();
norm_data = (image - min(image(:)))/(max(image(:))-min(image(:)));
set(gca,'YDir','normal');
imagesc(log10(abs(norm_data)),[-2.5,-1]); hold on % scaled display
%imagesc(log10(abs(norm_data))); hold on % non scaled
title(strcat('Nmodes: ',num2str(nModes_click),'- click on all positive copies first')); hold on 
    
% Get user input clicks
[x,y] = ginput(nModes_click*2); % Get user clicks for all modes
yp = y(1:nModes_click);
ym = y(nModes_click+1:end);
xp = x(1:nModes_click);
xm = x(nModes_click+1:end);
spotCoords_userInput = [yp,xp,ym,xm];
close();

% Refine user spot coords

spotCoords_px = refine_coords2(image, spotCoords_userInput, r_search, r_cent, centClip);



% Column subtract image (remove readout streak if present)
nModes_sense = size(spotCoords_px,1);

% Fit radial vectors to extracted regions
radCoeffs = zeros(2*nModes_sense,2); % Line coefficients (mode (+/- spots) , m/c)
intercept_Vectors = zeros((2*nModes_sense).^2,4); % Arracy for computing intercept estimate (all spot permutations)
%figure();
%subplot(1,2,1);
%imagesc(log10(image-min(image(:)))); hold on
%xVals = 1:size(image,2); hold on
%colorbar();
index = 1;
indices = [1 2*nModes_sense];
for i=1:nModes_sense
    
    % shorter wavelength
    [ROI_max,~,offset] = cropROI(image, spotCoords_px(i,1:2), r_search*4); % Extract ROI
    coeffs_plus = lstsqr_fit_img(ROI_max,offset); % Fit line to image data
    radCoeffs(index,:) = coeffs_plus;
    intercept_Vectors(indices(1):indices(2),3:end) = repmat(coeffs_plus,2*nModes_sense,1);
    index = index+1;
    indices = indices+2*nModes_sense;

    % longer wavelength
    [ROI_max,~,offset] = cropROI(image, spotCoords_px(i,3:4), r_search*4);
    coeffs_minus = lstsqr_fit_img(ROI_max,offset);
    radCoeffs(index,:) = coeffs_minus;
    intercept_Vectors(indices(1):indices(2),3:end) = repmat(coeffs_minus,2*nModes_sense,1);
    index = index+1;
    indices = indices+2*nModes_sense;
end

intercept_Vectors(:,1:2) = repmat(radCoeffs,nModes_sense*2,1); % Fill in first columns
    
% Calculate median intercept point x=(m2-m1)/(c1-c2), y=m1x+c1
x = (intercept_Vectors(:,4)-intercept_Vectors(:,2))./(intercept_Vectors(:,1)-intercept_Vectors(:,3));
y = intercept_Vectors(:,1).*x+intercept_Vectors(:,2);
dispersion_cent = [nanmedian(y),nanmedian(x)]; % Louis - can stop here if the estimates are good enough

dispersion_cent
plot(spotCoords_px(:,2),spotCoords_px(:,1),'bx'); hold on % Blue
plot(spotCoords_px(:,4),spotCoords_px(:,3),'rx'); hold on % red

Co_effs =  zeros(nModes_sense,2); % Line coefficients (mode (+/- spots) , m/c)


for i = 1:nModes_click
    xb = spotCoords_px(i,2);
    yb = spotCoords_px(i,1);
    xr = spotCoords_px(i,4);
    yr = spotCoords_px(i,3);
    m = (yb-yr)/(xb-xr);
    c  = yb-m*xb;
    Co_effs(i,:) = [m,c];
    x = 1:640;
    y = m*x+c;
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
dispersion_cent