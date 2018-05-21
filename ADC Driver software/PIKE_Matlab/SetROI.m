function  SetROI_maxval(vid,dim)

%re-set the ROI to full frame
vid.ROIPosition = [0 0 640 480];

% take exposure
image = take_exposure(vid);
% convert exposure to double 
image = im2double(image);

% get indecies of max val
[val,x] = max(max(image));
[val,y] = max(image(:,x));

% set boundary to dim x dim around center
ROIBoundary = [x-dim/2 y-dim/2 dim+1 dim+1];
vid.ROIPosition = ROIBoundary;


end