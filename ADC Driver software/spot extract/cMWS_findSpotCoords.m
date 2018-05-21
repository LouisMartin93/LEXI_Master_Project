% Finds coordinates of cMWS spots in pixel coordinates from scratch, using
% user clicks as input:
% - image: Calibration image for user alignment
% - spotCoords_lD: Full list of design spot coords in lambda/D (2,Nmodes,2)
% - nModes_click: Number of modes for the user will provide location
%                 (<= length of spotCoords_lD)
% NB - click on ALL I_p first, then ALL I_m spots.

function [spotCoords_px, transform, spotCoords_BB, zeroCoords_BB, dispersion_cent] = cMWS_findSpotCoords(image, spotCoords_lD, nModes_click, r_search, r_cent, centClip, bandwidth,lambda_cent)
    
    figure();
    norm_data = (image - min(image(:)))/(max(image(:))-min(image(:)));
    set(gca,'YDir','normal');
    %imagesc(log10(abs(norm_data)),[-3,-2]); hold on % scaled display
    imagesc(log10(abs(norm_data))); hold on % non scaled
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
    spotCoords_px = refine_coords(image, spotCoords_userInput, r_search, r_cent, centClip);
    % Find best fit to all spots based on user-supplied clicks
    [spotCoords_px,transform] = cMWS_updateSpotCoords(image,spotCoords_px,spotCoords_lD, r_search, r_cent, centClip);
    
    if bandwidth > 25
        for i =1:5
            % Iterate NB alignment, to ensure all coords fall somewhere along dispersed spots
            [spotCoords_px,transform] = cMWS_updateSpotCoords(image,spotCoords_px,spotCoords_lD, r_search, r_cent, centClip);
        end
        % Extract broadband spot coordinates (min,mean,max)
        [spotCoords_BB,zeroCoords_BB,dispersion_cent] = cMWS_updateSpotCoords_BB(image,spotCoords_px,spotCoords_lD,r_search,bandwidth,lambda_cent); 
    else
        spotCoords_BB = '';
        zeroCoords_BB = '';
        dispersion_cent = ''; 
    end
end

    
    
    