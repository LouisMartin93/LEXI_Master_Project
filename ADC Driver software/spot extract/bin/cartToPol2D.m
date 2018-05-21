% Cartesian to Polar coordinate conversion (2D) 
% NB: theta in radians 

function [r,theta] = cartToPol2D(y,x)
    
    r = sqrt(abs(y).^2 + abs(x).^2);
    theta = atan2(y,x);
    
end
    
