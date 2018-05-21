% Setup exopsure parameters
function dm = cont_exposures(N,video_obj,varargin)
    numvarargs = length(varargin);
    if numvarargs > 1
        error('cont_expousres.m: Too many outputs requires only 1 aditional output')
    end
    
    optargs = {'linear'};
    optargs(1:numvarargs) = varargin;
    DispMode = string(optargs);
    fig = figure('position',[200 200 400 400]);
    colormap 'gray'
    for i=1:N
        image = take_exposure(video_obj);
        image = im2double(image);
        if DispMode == 'log'
            imagesc(log10(image));
            title(max(image(:)));
            drawnow();
        end
        if DispMode == 'linear'
            imagesc(image);
            title(max(image(:)));
            drawnow();
        end
        

    end
end
