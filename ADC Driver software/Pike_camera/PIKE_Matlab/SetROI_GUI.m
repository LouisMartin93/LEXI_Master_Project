function  SetROI_GUI(vid,dim)

% take exposure
image = take_exposure(vid);
%display expousre only linear
image = im2double(image);
fig = figure;
imagesc(image);
% Get user input clicks
[x,y] = ginput(1); % Get user clicks for all modes
ROIBoundary = [x-dim/2 y-dim/2 dim dim];
vid.ROIPosition = ROIBoundary;
close(fig);

end