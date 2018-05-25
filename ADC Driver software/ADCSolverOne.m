function NewTheta = ADCSolverOne(Theta,P,dispersion_vector)
    %theta must be set correctly to the assumed zero point
    rx = dispersion_vector(2);
    ry = dispersion_vector(1);
    
    dThetax = real(acosd(rx/P) - Theta)
    dThetay = real(asind(ry/P) - Theta)
    dTheta = (dThetax +dThetay)/2
    
    NewTheta = Theta+dTheta;
    NewTheta = mod(NewTheta*180/pi,360);
    
end