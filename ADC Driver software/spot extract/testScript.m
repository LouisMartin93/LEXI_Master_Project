clear all;
close all;

addpath('bin');

% Make simulated PSF
%aper = fitsread('cMWS_360\aper.fits');
%APP_phase = fitsread('cMWS_360\cMWS.fits');
%NCPE_phase = zeros(size(APP_phase));
lD2px = 3.2; % pixels per lambda/D
imSize = [2*640,2*480]; % pixels
bandwidth = 145; %nm
lambda_cent = 560.5; % Central wavelength (nm)

%image = make_simPSF(aper, APP_phase, NCPE_phase, lD2px, imSize, bandwidth, lambda_cent);


%if using own data
image = fitsread('../PIKE_Matlab/Test Images/dispersion/SM450/f4_cMWS/image1.fits');

imagesc(log10(image)); hold on
colorbar(); hold on


% Find spots
spotCoords_lD = fitsread('cMWS_360\spotCoords.fits');
nModes_click = 1; % Number of modes to manually select
r_search = 5*lD2px; % ROI for finding spot maxima (px)
r_cent = 3*lD2px; % ROI for calculating centroid (px)
centClip = 5; % Lower clip threshold for centroiding (peak flux / centClip)

[spotCoords_px, transform, spotCoords_BB, zeroCoords_BB, dispersion_cent] = cMWS_findSpotCoords(image, spotCoords_lD, nModes_click, r_search, r_cent, centClip, bandwidth, lambda_cent);

% Plot coordinates
if bandwidth >25
        % Spot coordinates
        plot(spotCoords_BB(:,2,1),spotCoords_BB(:,1,1),'bx'); hold on % Blue, +
        plot(spotCoords_BB(:,4,1),spotCoords_BB(:,3,1),'bx'); hold on % Blue, -
        plot(spotCoords_BB(:,2,2),spotCoords_BB(:,1,2),'gx'); hold on % Green, +
        plot(spotCoords_BB(:,4,2),spotCoords_BB(:,3,2),'gx'); hold on % Green, -
        plot(spotCoords_BB(:,2,3),spotCoords_BB(:,1,3),'rx'); hold on % Red, +
        plot(spotCoords_BB(:,4,3),spotCoords_BB(:,3,3),'rx'); hold on % Red, -
        % Central coordinates
        plot(zeroCoords_BB(2,1),zeroCoords_BB(1,1),'bx'); hold on % Blue, -
        plot(zeroCoords_BB(2,2),zeroCoords_BB(1,2),'gx'); hold on % Green, -
        plot(zeroCoords_BB(2,3),zeroCoords_BB(1,3),'rx'); hold on % Red, -
        % Dispersion center
        plot(dispersion_cent(2),dispersion_cent(1),'kx');
else
        plot(spotCoords_px(:,2),spotCoords_px(:,1),'bx'); hold on % +
        plot(spotCoords_px(:,4),spotCoords_px(:,3),'bx'); hold on % -
        % Central coords
        plot(transform(2),transform(1),'kx');
end
