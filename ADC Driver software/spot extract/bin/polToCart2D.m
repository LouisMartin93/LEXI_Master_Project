% Polar to Cartesian coordinate conversion (2D) 
% NB: theta in radians 

function [y,x] = polToCart2D(r,theta)
    
    y = r.*sin(theta);
    x = r.*cos(theta);
    
end