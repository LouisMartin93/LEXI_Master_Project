% ADC calibration for a single prism ADC
% based on the works by Pathek et al.

function P = ADCCalibOne(ADCHandle,Theta,ThetaOffset,videoObject,DarkFrame)
l = cosd(Theta)  - cosd(Theta-ThetaOffset);
m = sind(Theta)  - sind(Theta-ThetaOffset);

%move ADC to positions
ADCMovePos(ADCHandle,Theta)
pause(2); % wait for motors to be in positions
%next comes the image aquistion and check that it is aligned correctly as
%well as calculates the dispersion
dispersion_vector1 = 2*OneIteration(videoObject,DarkFrame);

ADCMovePos(ADCHandle,Theta-ThetaOffset);
pause(2);
dispersion_vector2 = 2*OneIteration(videoObject,DarkFrame);

%calculates the prisms response P
P = ((dispersion_vector1(2)-dispersion_vector2(2))^2 + (dispersion_vector1(1)-dispersion_vector2(1))^2)/(l^2 +m^2);

end