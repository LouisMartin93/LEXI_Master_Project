% Calculates the transform between lambda/D and pixel coordinates from the
% given arrays (do not have to be of equal length, but are of the form (mode,y/x/p/m)
% (i.e. of size (Nmodes,4)) - both in x/y coordinates. rotAngle is 
% computed clockwise, from lD to px.

function [newCoords_px, transform] = spotCoords_lD2px(spotCoords_px_0,spotCoords_lD_0)

    % Re-format coords to be of form (Nmodes*2,2)
    Nmodes_px = size(spotCoords_px_0,1);
    Nmodes_lD = size(spotCoords_lD_0,1);
    spotCoords_px = zeros(Nmodes_px*2,2);
    spotCoords_lD = zeros(Nmodes_lD*2,2);
    spotCoords_px(1:Nmodes_px,:) = spotCoords_px_0(:,1:2);
    spotCoords_px(Nmodes_px+1:end,:) = spotCoords_px_0(:,3:4);
    spotCoords_lD(1:Nmodes_lD,:) = spotCoords_lD_0(:,1:2);
    spotCoords_lD(Nmodes_lD+1:end,:) = spotCoords_lD_0(:,3:4);

    % Compute PSF center
    PSF_cent = nanmean(spotCoords_px,1);
    %plot(PSF_cent(2),PSF_cent(1),'wx')
    
    [rs_lD,thetas_lD] = cartToPol2D(spotCoords_lD(:,1),spotCoords_lD(:,2));
    [rs_px,thetas_px] = cartToPol2D(spotCoords_px(:,1)-PSF_cent(1),spotCoords_px(:,2)-PSF_cent(2));
    
    % Correct thetas with opposite sign (phase wrapping)
    %ind_p = find(thetas_lD(1:Nmodes_px)>pi/2 & thetas_px(1:Nmodes_px)<-pi/2); % Positive spots
    %ind_m = find(thetas_px(1:Nmodes_px)>pi/2 & thetas_lD(1:Nmodes_px)<-pi/2);
    %thetas_px(ind_p) =  2*pi - thetas_px(ind_p);
    %thetas_lD(ind_m) =  2*pi - thetas_lD(ind_m);
    %ind_p = find(thetas_lD(Nmodes_lD+1:Nmodes_lD+Nmodes_px)>pi/2 & thetas_px(Nmodes_px+1:end)<-pi/2); % Negative spots
    %ind_m = find(thetas_px(Nmodes_px+1:end)>pi/2 & thetas_lD(Nmodes_lD+1:Nmodes_lD+Nmodes_px)<-pi/2);
    %thetas_px(ind_p) =  2*pi - thetas_px(ind_p);
    %thetas_lD(ind_m) =  2*pi - thetas_lD(ind_m);
    
    % Compute mean radius conversion (i.e. lambda/D to px) for available points
    rotAngle = [(thetas_px(1:Nmodes_px) - thetas_lD(1:Nmodes_px)) (thetas_px(Nmodes_px+1:end) - thetas_lD(Nmodes_lD+1:Nmodes_lD+Nmodes_px))];
    rotAngle(rotAngle>pi) = rotAngle(rotAngle>pi)-2*pi;
    rotAngle(rotAngle<pi) = rotAngle(rotAngle<pi)+2*pi;
    rotAngle(rotAngle>2*pi) = rotAngle(rotAngle>2*pi)-2*pi;
    rotAngle(rotAngle<-2*pi) = rotAngle(rotAngle<-2*pi)+2*pi;
    rotAngle = nanmedian(rotAngle(:)); % Rotation angle
   
    % Compute mean lambda/D -> pixel scaling for available points
    lD2px = [(rs_px(1:Nmodes_px)./rs_lD(1:Nmodes_px)) (rs_px(Nmodes_px+1:end)./rs_lD(Nmodes_lD+1:Nmodes_lD+Nmodes_px))];
    lD2px = nanmean(lD2px(:));

    % Convert lD coordinates to px
    thetas_lD2px = thetas_lD+rotAngle;
    rs_lD2px = rs_lD.*lD2px;
    [y_lD2px,x_lD2px] = polToCart2D(rs_lD2px,thetas_lD2px);
    y_lD2px = y_lD2px+PSF_cent(1);
    x_lD2px = x_lD2px+PSF_cent(2);
    newCoords_px = zeros(Nmodes_lD,4);
    newCoords_px(:,1) = y_lD2px(1:Nmodes_lD);
    newCoords_px(:,2) = x_lD2px(1:Nmodes_lD);
    newCoords_px(:,3) = y_lD2px(Nmodes_lD+1:end);
    newCoords_px(:,4) = x_lD2px(Nmodes_lD+1:end);
    
    transform = [PSF_cent(1),PSF_cent(2),lD2px,rotAngle];
    
end
