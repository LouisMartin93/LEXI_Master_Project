% Returns (smoothed) maximum and centroid WFS spot coordinates for supplied image
% NB - assumes coords are of the form [[yp,xp,ym,xm],...]

function  [newCoords] = refine_coords(image, oldCoords, r_search, r_cent, centClip)

    newCoords = zeros(size(oldCoords));

    for i = 1:length(oldCoords(:,1))
        
        % I_plus
        [ROI_max,centCoords_zoom] = cropROI(image, oldCoords(i,1:2), r_search*2);
        coords_max = maxInd(ROI_max)+ oldCoords(i,1:2)-centCoords_zoom; % Max 
        [ROI_cent,centCoords_zoom] = cropROI(image, coords_max, r_cent*2);
        coords_plus = centroid(ROI_cent,centClip); % Centroid
        
        %subplot(2,2,1)
        %imagesc(ROI_cent); hold on
        %plot(coords_plus(2),coords_plus(1),'kx')
        
        coords_plus = coords_plus + coords_max-centCoords_zoom;
        newCoords(i,1:2) = coords_plus;
       

        %subplot(2,2,2)
        %imagesc(log10(abs(image))); hold on
        %plot(coords_plus(2),coords_plus(1),'kx')

        % I_minus
        [ROI_max,centCoords_zoom] = cropROI(image, oldCoords(i,3:4), r_search*2);
        coords_max = maxInd(ROI_max)+ oldCoords(i,3:4)-centCoords_zoom; % Max 
        [ROI_cent,centCoords_zoom] = cropROI(image, coords_max, r_cent*2);
        coords_minus = centroid(ROI_cent,centClip); % Centroid
        
        %subplot(2,2,3)
        %imagesc(ROI_cent); hold on
        %plot(coords_minus(2),coords_minus(1),'kx')
        
        coords_minus = coords_minus+ coords_max-centCoords_zoom;
        newCoords(i,3:4) = coords_minus;
        
        %subplot(2,2,4)
        %imagesc(log10(abs(image))); hold on
        %plot(coords_minus(2),coords_minus(1),'kx')
        %pause()

        if sum(isnan(newCoords(i,:)))>0
            % Delete entire mode coordinates
            newCoords(i,:) = nan;
        end
            
    %close()
    end
end
        
