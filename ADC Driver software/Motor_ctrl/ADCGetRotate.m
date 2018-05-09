%% Returns the rotation angle (position) need in the ADC for a given zenith
function [rotation_position] = ADCGetRotate(zenith)
    a=0.5085; % calibraion factor previously calculated
    rotation_position = acos(-a*tan(zenith*pi/180))*180/pi - 90;% rotation needed
end