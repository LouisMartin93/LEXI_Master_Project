function disp_exposure(image,varargin)
    % Displays an image of the exposure must first convert from 
    % uint8 to a double.
    numvarargs = length(varargin);
    if numvarargs > 2
        error('disp_expousres.m: Too many outputs requires only 2 aditional output')
    end
    
    optargs = {'both',true};
    optargs(1:numvarargs) = varargin;
    [DispMode,norm] = optargs{:};
    
    
    %normalising image 
    image = im2double(image);
    if norm
        image = image./max(image(:));
    end
    % displaying subplots of linear and log10 format of image for
    % comparison.
    if string(DispMode)== 'both'
        figure('Position', [125,300,1000,330]);
        subplot(121);
        hold on;
        imagesc(log10(image));
        title('log10 scale')
        colormap gray;
        colorbar();
        axis tight
        hold off
        subplot(122);
        hold on;
        imagesc(image);
        title('linear scale')
        colormap gray;
        colorbar();
        axis tight
        hold off
    end
    if string(DispMode)== 'linear'
        figure('Position', [125,300,500,400]);
        imagesc(image);
        title('linear scale')
        colormap gray;
        colorbar();
        axis tight
        hold off
    end 
    if string(DispMode)== 'log'
        figure('Position', [125,300,500,480]);
        imagesc(log10(image));
        title('linear scale')
        colormap gray;
        colorbar();
        axis tight
        hold off
    end
end
