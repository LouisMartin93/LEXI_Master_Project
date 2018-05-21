% Finds coordinates of maximum pixel of an NxN pixel region

function [maxCoords,maxVal] = maxInd(region)

    [maxVals,y] = max(region);
	[maxVal,x] = max(maxVals);
	maxCoords = [y(x),x];

end