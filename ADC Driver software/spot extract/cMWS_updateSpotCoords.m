% Converts theoretical cMWS design spot coords (in lambda/D) to pixels
% based on the scaling, rotation fit of (a subset of) pixel coords
% - spotCoords_px: List of pixel coordinates (may be a subset of total
% - spotCoords_lD: The full list of design coordinates
% - transform: [centY,centX,scaleY,scaleY,rotAngle] (pixels,ratio,radians)

function [spotCoords_px, transform] = cMWS_updateSpotCoords(image,spotCoords_px_0, spotCoords_lD, r_search, r_cent,centClip)

    spotCoords_px = spotCoords_px_0;

    if size(spotCoords_px_0,1) ~= size(spotCoords_lD,1)
        % If an incomplete subset of spotCoords_lD is supplied in
        % spotCoords_px, extend to predict all points
        [spotCoords_px, ~] = spotCoords_lD2px(spotCoords_px,spotCoords_lD);   
    end
    
    % Re-compute centroid
    spotCoords_px = refine_coords(image, spotCoords_px, r_search, r_cent,centClip);
    
    % Refine template alignment
    [spotCoords_px, transform] = spotCoords_lD2px(spotCoords_px,spotCoords_lD);
        
end