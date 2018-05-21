% Crops the specified image to the desired dimensions, returning the ROI

function [ROI,centCoords_zoom,crop_offset] = cropROI(img, cropCent, cropSize)
    
    num_x_pixels = size(img,1);
    num_y_pixels = size(img,2);

    if length(cropSize) == 1
        cropSize = [cropSize,cropSize];
    end
    cropCent(isnan(cropCent)) = -cropSize(1)*2;
        
    cropCoords = [cropCent(1)-cropSize(1)/2,cropCent(1)-cropSize(1)/2+cropSize(1)-1, ...
                  cropCent(2)-cropSize(2)/2,cropCent(2)-cropSize(2)/2+cropSize(2)-1];
              
    
	if cropCoords(1) > num_y_pixels || cropCoords(3) > num_x_pixels || cropCoords(2) < 1 || cropCoords(4) < 1
        % Check for completely incorrect ROIs outside FoV
        %disp('Warning: Tried to crop a region completely outside the image! Returning array of nans') 
        ROI = nan(num_y_pixels,num_x_pixels);
        centCoords_zoom = [nan,nan];
        crop_offset = [nan,nan];
    else
        % Check for indices outside array
        if cropCoords(1) < 1
            cropCoords(1) = 1;
        end
        if cropCoords(3) < 1
            cropCoords(3) = 1;
        end
        if cropCoords(2) > num_y_pixels
            cropCoords(2) = num_y_pixels;
        end
        if cropCoords(4) > num_x_pixels
            cropCoords(4) = num_x_pixels;
        end

        %cropSize
        %cropCoords
        %cropCoords(2)-cropCoords(1)
        %cropCoords(4)-cropCoords(3)

        cropCoords_round = int32(round(cropCoords));   
        ROI = img(cropCoords_round(1):cropCoords_round(2),cropCoords_round(3):cropCoords_round(4));

        % Calculate (to avoid rounding error)
        crop_offset = double([cropCoords_round(1),cropCoords_round(3)]);
        centCoords_zoom = cropCent - crop_offset +1;
    end


end
