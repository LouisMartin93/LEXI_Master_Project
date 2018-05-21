function [ image ] = make_simPSF(aper, APP_phase, NCPE_phase, lD2px, imSize, bandwidth, lambda_cent)
    % Generates broadband, noisy, dispersed vAPP PSFs of the supplied phase pattern
    % aper (NxN) - vAPP pupil mask
    % APP_phase (NxN) - vAPP phase 
    % NCPE_phase (NxN) - aberrating phase
    % lD2px (float) - image sampling (pixel/lD)
    % imSize (float) - image size (pixels)
    % bandwidth (float) - image bandwidth (nm)

    phase_plus = APP_phase + NCPE_phase;
	phase_minus = -APP_phase + NCPE_phase;
	phase_zero = NCPE_phase;
	padSize = round(size(aper,1).*(lD2px-1)./2.);
	num_x_pixels = imSize(2);
	num_y_pixels = imSize(1);
    
    circPolFrac = 1;%0.2; % Fractional excess of one circular polarisation
    leakage = 5./100.; % Fractional vAPP leakage
    
    % Make PSF (monochromatic) then if necessary re-sample for BB.
	wf_plus=aper.*exp(1.j.*phase_plus);
	wf_minus=aper.*exp(1.j.*phase_minus);
	wf_zero=aper.*exp(1.j.*phase_zero);
    PSF_0 = abs(FT(wf_plus,padSize)).^2+(abs(FT(wf_minus,padSize)).^2)*circPolFrac+(abs(FT(wf_zero,padSize)).^2)*leakage;
    if bandwidth == 0
        [image,~] = cropROI(PSF_0, round(size(PSF_0)./2), imSize);
    else
        Nsamples = 31; % Number of discrete wavelengths to simulate
        disp_ldnm = 0.02; % Dispersion rate (lambda/D per nm)

        disp_pxnm = disp_ldnm*lD2px; % Dispersion (pixels/nm)
        lambdas = linspace(lambda_cent-bandwidth/2.,lambda_cent+bandwidth/2.,Nsamples);
        sampling = lD2px*lambdas/lambda_cent; % Sampling (in px/lD) for each wavelength
        PSF_shifts = int32((lambdas-lambda_cent).*disp_pxnm); % Pixel shifts due to dispersion
        image = zeros(num_y_pixels,num_x_pixels);
        for i = 1:Nsamples
            % Scale and crop image (may have 1px shift error)
            image_i = imresize(PSF_0,sampling(i)/lD2px);
            cropCoords = int32([(size(image_i,1)-num_y_pixels)/2,(size(image_i,2)-num_x_pixels)/2.]);
            image_i = image_i(cropCoords(1):cropCoords(1)+num_y_pixels-1,cropCoords(2):cropCoords(2)+num_x_pixels-1);
            image_i = circshift(image_i,PSF_shifts(i),2);
            image = image+image_i;
        end
    end
    
    % Add noise and readout line
    noiseContrast = 5.E-4; % RO-noise rms wrt. peak flux
    peakFlux = 12000; % Counts
    bias = 400; % Bias offset
    ROline_contrast = 1E-3; % Contrast of readout striping line
    
    colSum = max(image);
	image = image + colSum*ROline_contrast; % Readout line
    noise = normrnd(0,sqrt(image)); % Photon noise
    noise = noise + normrnd(0,max(image(:)).*noiseContrast,size(image,1),size(image,2)); % RO noise
	image = image+noise;
    image = image./max(image(:)).*peakFlux+bias; % Convert to counts + bias
	image = double(int32(image));
    
    image = flipud(image);
    image = image.';

    
end