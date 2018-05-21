% Least-squared fitting for line in image data

function [coeffs] = lstsqr_fit_img(image,offset)
    
    cutoff = 0.5; % Contrast cut
    ii = find(image>max(image(:))*cutoff);
    [y,x] = ind2sub(size(image),ii);
    y = y+offset(1);
    x = x+offset(2);
    [coeffs,S] = polyfit(x,y,1);
    [ys,d] = polyval(coeffs,x,S);
    
    % Compute flipped line of flipped structure (polyfit doesn't work
    % with vertical lines)
    [coeffs_flip,S_flip] = polyfit(y,x,1);
    %coeffs_flip = flipline(coeffs_flip(1),coeffs_flip(2),size(image,1)/2);
    [ys_flip,d_flip] = polyval(coeffs_flip,y,S_flip);
    
    rms = std(d);
    rms_flip = std(d_flip);
    if rms_flip<rms
        % Replace with final plot
        coeffs(1) = 1/coeffs_flip(1);
        coeffs(2) = -coeffs_flip(2)/coeffs_flip(1);
    end
    
