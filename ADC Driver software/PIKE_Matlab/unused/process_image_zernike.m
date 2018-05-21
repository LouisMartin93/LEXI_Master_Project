function [phase] = process_image_zernike(data)
    
    I_left = select_pupil_data(data, center, radius);
    I_right = select_pupil_data(data, center, radius);
    
    phase = (I_left-I_right)./(I_left + I_right);
 
end

