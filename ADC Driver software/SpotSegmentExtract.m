%% Function for New image processing technique for dispersion centre extracion
%
%
% This is built up from the script quickScript_imgProcess.m which follows
% similar guidelines as done in the Otsu method, see https://uk.mathworks.com/help/images/examples/correcting-nonuniform-illumination.html
%
% Input 
% a dark subtracted image, that is all that needs to be done and the
%
% Output
% the Dispersion vector.

%% Begin Function

function [dispersion_vector,dispersion_cent,centre] = SpotSegmentExtract(image,varargin)

    %% Checking and seting Varargin
    numvarargs = length(varargin);
    if numvarargs > 1
        error('SpotSegmentExtract.m: Too many outputs requires only 1 aditional output')
    end
    
    optargs = {false};
    optargs(1:numvarargs) = varargin;
    PlotPlots = optargs{:};

    
    %% normalise image
    image = image - min(image(:));
    
    %% gaussian filter image
    I2 = imgaussfilt(image,3);
    if PlotPlots == true
        figure;imagesc(I2);
        title('sigma = 3')
    end
    %% enhance contrast
    I3 = imadjust(I2);
    if PlotPlots == true
        figure(); imagesc(I3);
        title('sigma = 3')
    end
    
    %% binarize image 
    bw = imbinarize(I3);
    if PlotPlots == true
        figure;
        imagesc(bw);
    end
    %% get rid of backround binary noise
    bw1 = bwareaopen(bw,100);
    if PlotPlots == true
        figure;imagesc(bw1);
    end
    %% Set the spot properties table up and remove central region.
    cc = bwconncomp(bw1,4);
    spotsdata = regionprops('table',cc,'Area','Centroid','Eccentricity','Orientation');

    spotareas = [spotsdata.Area];
    [max_area,idx] = max(spotareas);
    centre1 = [spotsdata.Centroid(idx,2) spotsdata.Centroid(idx,1)];
    spotsdata(idx,:) =[];

    %% removing the small spots which are likely to throw the centering.
    spotareas = [spotsdata.Area];
    med = median(spotareas);
    standev = std(spotareas);
    idxs = find(spotareas < med-2*standev);
    spotsdata(idxs,:) =[];

    %% making plots without the central spot (central PSF)
    spots = false(size(bw1));
    for i = 1:cc.NumObjects
        spots(cc.PixelIdxList{i}) = true;
    end
    spots(cc.PixelIdxList{idx}) = false;
    
    if PlotPlots == true
        figure;
        imagesc(spots);
        hold on
    end
    
    %% Calculating the intercept and slope of each Orientation of each spot.
    x = 1:640;
    Co_effs =  zeros(size(spotsdata,1),2);
    for i = 1:size(spotsdata,1)
        theta = spotsdata.Orientation(i);
        x1 = spotsdata.Centroid(i,1);
        y1 = spotsdata.Centroid(i,2);
        m = -tand(theta);
        c = y1+tand(theta)*x1;
        Co_effs(i,:) = [m,c];
        if PlotPlots == true
            plot(x,-x*tand(theta) + c,'-r')
            hold on
        end
    end
    
    %% Calculating the Median of the intercepts of each line vectors for the dispersion centre.
    
    numSpots = size(spotsdata,1);
    xs = zeros(1,3);
    ys = zeros(1,3);
    c = combnk(1:numSpots,2);
    for i= 1:size(c,1)
        ccell = num2cell(c(i,:));
        [id1,id2] = deal(ccell{:});
        xs(i) = (Co_effs(id2,2)- Co_effs(id1,2))/(Co_effs(id1,1)- Co_effs(id2,1));
        ys(i) = Co_effs(id2,1)*xs(i) + Co_effs(id2,2);
    end
    
    %% calculating the PSF center
    % This isnt too far off the spot centre of the PSF (1-2 pixels).
    centre = Centroid2DGaussian(image);
    
    dispersion_cent = [nanmedian(ys),nanmedian(xs)];
    dispersion_vector = dispersion_cent - centre;
    
    if PlotPlots == true
        plot(dispersion_cent(2),dispersion_cent(1),'kx');
        hold off
    
        figure();
        imagesc(log10(abs(image)));
        hold on
        plot(dispersion_cent(2),dispersion_cent(1),'kx');
        plot(centre(2),centre(1),'k+')
    end 
end
