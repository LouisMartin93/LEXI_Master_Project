% Identifies broadband spot extraction coordinates (min, mean, max) from
% input image and full set of centroided spot coords (pxels)
% - Identifies radial vector for each spot and computes mean intersection
% - Makes SNR cut along new radial vectors and original points to identify
%   edges and central wavelength bin
% - Uses dispersion model and lD separations to refine min/max positions

function [spotCoords_BB,zeroCoords_BB,dispersion_cent] = cMWS_updateSpotCoords_BB(image,spotCoords_px,spotCoords_lD,r_search,bandwidth,lambda_cent)

    % Column subtract image (remove readout streak if present)
    image = image - median(image(1:50,:));
    nModes_sense = size(spotCoords_lD,1);

    % Fit radial vectors to extracted regions
    radCoeffs = zeros(2*nModes_sense,2); % Line coefficients (mode (+/- spots) , m/c)
    intercept_Vectors = zeros((2*nModes_sense).^2,4); % Array for computing intercept estimate (all spot permutations)
    %figure();
    %subplot(1,2,1);
    %imagesc(log10(image-min(image(:)))); hold on
	%xVals = 1:size(image,2); hold on
	%colorbar();
    index = 1;
    indices = [1 2*nModes_sense];
    for i=1:nModes_sense
        % I_plus spot
        [ROI_max,~,offset] = cropROI(image, spotCoords_px(i,1:2), r_search*4); % Extract ROI
        coeffs_plus = lstsqr_fit_img(ROI_max,offset); % Fit line to image data
        radCoeffs(index,:) = coeffs_plus;
        intercept_Vectors(indices(1):indices(2),3:end) = repmat(coeffs_plus,2*nModes_sense,1);
        index = index+1;
        indices = indices+2*nModes_sense;

        % I_minus spot
        [ROI_max,~,offset] = cropROI(image, spotCoords_px(i,3:4), r_search*4);
        coeffs_minus = lstsqr_fit_img(ROI_max,offset);
        radCoeffs(index,:) = coeffs_minus;
        intercept_Vectors(indices(1):indices(2),3:end) = repmat(coeffs_minus,2*nModes_sense,1);
        index = index+1;
        indices = indices+2*nModes_sense;
    end

    intercept_Vectors(:,1:2) = repmat(radCoeffs,nModes_sense*2,1); % Fill in first columns
    
    % Calculate median intercept point x=(m2-m1)/(c1-c2), y=m1x+c1
	x = (intercept_Vectors(:,4)-intercept_Vectors(:,2))./(intercept_Vectors(:,1)-intercept_Vectors(:,3));
    y = intercept_Vectors(:,1).*x+intercept_Vectors(:,2);
    dispersion_cent = [nanmedian(y),nanmedian(x)]; % Louis - can stop here if the estimates are good enough
    
    % Re-compute radial vectors for each originally detected spot coord
    radCoeffs = zeros(size(spotCoords_px));
    radCoeffs(:,1) = (spotCoords_px(:,1)-dispersion_cent(1))./(spotCoords_px(:,2)-dispersion_cent(2));
    radCoeffs(:,2) = dispersion_cent(1)-radCoeffs(:,1).*dispersion_cent(2);
    radCoeffs(:,3) = (spotCoords_px(:,3)-dispersion_cent(1))./(spotCoords_px(:,4)-dispersion_cent(2));
    radCoeffs(:,4) = dispersion_cent(1)-radCoeffs(:,3).*dispersion_cent(2);
    
    %for i =1:nModes_sense
    %    plot(xVals,xVals*radCoeffs(i,1)+radCoeffs(i,2),'w'); hold on
    %    plot(xVals,xVals*radCoeffs(i,3)+radCoeffs(i,4),'w'); hold on
    %    plot(spotCoords_px(:,2),spotCoords_px(:,1),'rx'); hold on
    %    plot(spotCoords_px(:,4),spotCoords_px(:,3),'rx'); hold on
    %end
    
    %plot(dispersion_cent(2),dispersion_cent(1),'rx');
    
    % Identify spot minima/maxima
    rs = zeros(size(spotCoords_px,1),2); % Radial distances of spots from dispersion center
    thetas = zeros(size(spotCoords_px,1),2);
    rs(:,1) = sqrt(abs(spotCoords_px(:,1)-dispersion_cent(1)).^2+abs(spotCoords_px(:,2)-dispersion_cent(2)).^2);
    rs(:,2) = sqrt(abs(spotCoords_px(:,3)-dispersion_cent(1)).^2+abs(spotCoords_px(:,4)-dispersion_cent(2)).^2);
    thetas(:,1) = atan2(spotCoords_px(:,1)-dispersion_cent(1),spotCoords_px(:,2)-dispersion_cent(2));
    thetas(:,2) = atan2(spotCoords_px(:,3)-dispersion_cent(1),spotCoords_px(:,4)-dispersion_cent(2));
    % Calculate dispersion size estimates given radial separation
    seg_lengths = rs*bandwidth/lambda_cent; % Fractional bandwidth = dispersion

    %figure();
    spotCoords_BB = zeros(nModes_sense,4,3);
    for index = 1:nModes_sense
        
        % Positive spot
        xs = [spotCoords_px(index,2)-abs(seg_lengths(index,1)*cos(thetas(index,1))) spotCoords_px(index,2)+abs(seg_lengths(index,1)*cos(thetas(index,1)))];
        ys = [radCoeffs(index,1)*xs(1)+radCoeffs(index,2) radCoeffs(index,1)*xs(2)+radCoeffs(index,2)];
        [cx,cy,c] = improfile(image,xs,ys);
        c(isnan(c)) = 0;
        c = medfilt1(c,5);
        edges = findchangepts(c,'MaxNumChanges',2); % Edge detection
        if isempty(edges)
            edges_y = [nan,nan];
            edges_x = [nan,nan];
        else
            edges_y = cy(edges).'.*ones(1,2);
            edges_x = cx(edges).'.*ones(1,2);
        end
        if max(edges_x) < dispersion_cent(2)
            edges_x = fliplr(edges_x);
            edges_y = fliplr(edges_y);
        end
        spotCoords_BB(index,1:2,1) = [edges_y(1),edges_x(1)];
        spotCoords_BB(index,1:2,3) = [edges_y(2),edges_x(2)];
        spotCoords_BB(index,1:2,2) = [mean(edges_y) mean(edges_x)];
        
        % Negative spot
        xs = [spotCoords_px(index,4)-abs(seg_lengths(index,2)*cos(thetas(index,2))) spotCoords_px(index,4)+abs(seg_lengths(index,2)*cos(thetas(index,2)))];
        ys = [radCoeffs(index,3)*xs(1)+radCoeffs(index,4) radCoeffs(index,3)*xs(2)+radCoeffs(index,4)];
        [cx,cy,c] = improfile(image,xs,ys);
        c(isnan(c)) = 0;
        c = medfilt1(c,5);
        edges = findchangepts(c,'MaxNumChanges',2);
        if isempty(edges)
            edges_y = [nan,nan];
            edges_x = [nan,nan];
        else
            edges_y = cy(edges).'.*ones(1,2);
            edges_x = cx(edges).'.*ones(1,2);
        end
        %edges_y = [cy(edges(1)) cy(edges(2))];
        %edges_x = [cx(edges(1)) cx(edges(2))];
        if max(edges_x) < dispersion_cent(2)
            edges_x = fliplr(edges_x);
            edges_y = fliplr(edges_y);
        end
        spotCoords_BB(index,3:4,1) = [edges_y(1),edges_x(1)];
        spotCoords_BB(index,3:4,3) = [edges_y(2),edges_x(2)];
        spotCoords_BB(index,3:4,2) = [mean(edges_y) mean(edges_x)];
        
        if sum(isnan(spotCoords_BB(index,:,1))) >0
            spotCoords_BB(index,:,1) = nan;
        end
        if sum(isnan(spotCoords_BB(index,:,2))) >0
            spotCoords_BB(index,:,2) = nan;
        end
        if sum(isnan(spotCoords_BB(index,:,3))) >0
            spotCoords_BB(index,:,3) = nan;
        end
        
        %clf();
        %plot(c); hold on
        %plot(spotCoords_px(index,2)-min(xs),0,'rx');
        %plot([edges(1),edges(1)],[0,max(c)],'k')
        %plot([edges(2),edges(2)],[0,max(c)],'k')
        %plot([centPx,centPx],[0,max(c)],'k')
        %pause();
    end
    
    % Separately fit RGB coeffs to model
    [spotCoords_BB(:,:,1), ~] = spotCoords_lD2px(spotCoords_BB(:,:,1),spotCoords_lD);
    [spotCoords_BB(:,:,2), ~] = spotCoords_lD2px(spotCoords_BB(:,:,2),spotCoords_lD);
    [spotCoords_BB(:,:,3), ~] = spotCoords_lD2px(spotCoords_BB(:,:,3),spotCoords_lD);
    
    %plot(spotCoords_BB(:,2,3),spotCoords_BB(:,1,3),'rx');hold on
    %plot(spotCoords_BB(:,4,3),spotCoords_BB(:,3,3),'rx');hold on
    %plot(spotCoords_BB(:,2,2),spotCoords_BB(:,1,2),'gx');hold on
    %plot(spotCoords_BB(:,4,2),spotCoords_BB(:,3,2),'gx');hold on
    %plot(spotCoords_BB(:,2,1),spotCoords_BB(:,1,1),'bx');hold on
    %plot(spotCoords_BB(:,4,1),spotCoords_BB(:,3,1),'bx');hold on
    
    % Compute PSF center for each RGB wavelength
    zeroCoords_BB = zeros(2,3);
    zeroCoords_BB(1,1) = mean(cat(1,spotCoords_BB(:,3,1),spotCoords_BB(:,1,1)));
    zeroCoords_BB(2,1) = mean(cat(1,spotCoords_BB(:,4,1),spotCoords_BB(:,2,1)));
    zeroCoords_BB(1,2) = mean(cat(1,spotCoords_BB(:,3,2),spotCoords_BB(:,1,2)));
    zeroCoords_BB(2,2) = mean(cat(1,spotCoords_BB(:,4,2),spotCoords_BB(:,2,2)));
    zeroCoords_BB(1,3) = mean(cat(1,spotCoords_BB(:,3,3),spotCoords_BB(:,1,3)));
    zeroCoords_BB(2,3) = mean(cat(1,spotCoords_BB(:,4,3),spotCoords_BB(:,2,3)));
    
    %plot(zeroCoords_BB(2,1),zeroCoords_BB(1,1),'bx');hold on
    %plot(zeroCoords_BB(2,2),zeroCoords_BB(1,2),'gx');hold on
    %plot(zeroCoords_BB(2,3),zeroCoords_BB(1,3),'rx');hold on
        
end
