% Centroids NxN pixel region around stated coordinate center, returning 
% sub-pixel-accurate coordinates

function [centCoords] = centroid(region,clip)
    
    [X,Y] = meshgrid(1:length(region(1,:)),1:length(region(:,1)));

    %arrayMed = median(median(region));
    arrayMax = max(max(region));
    region_sub = region - arrayMax./clip; %-arrayMed% Contrast of <=1/clip contributes to centroid
    region_sub(region_sub <= 0) = 0; % Clip to >= zero
    %figure(3)
    %imagesc(region_sub);hold on 
    %colorbar()       
    
    arraySum = sum(sum(region_sub));
    xWeight = sum(sum(region_sub.*X));
    yWeight = sum(sum(region_sub.*Y));
    
    centCoords = [yWeight/arraySum,xWeight/arraySum];
    %plot(centCoords(2),centCoords(1),'wx');
    %pause(0.5)
    %close()

end