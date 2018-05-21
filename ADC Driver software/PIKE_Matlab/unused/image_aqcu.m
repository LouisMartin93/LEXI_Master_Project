%clear work space
clear all;
close all;


% initialise camera that will be used.
vid_object = init_camera();

%source Properites for the exposure;
ROI = [314 148 100 100];
Shutter = 1.0;
Brightness = 128;
AutoExpo = 125;
Gain =  0;

% set source Properties
src = getselectedsource(vid_object);
vid_object.ROIPosition = ROI;
src.Shutter = Shutter;
src.Brightness = Brightness;
src.AutoExposure = AutoExpo;

% taking an exposure
start_camera(vid_object);
image = take_exposure(vid_object);
stop_camera(vid_object);

%normalising image 
image = im2double(image);
image = image./max(image(:));

% displaying log10 format
figure;
hold on;
imagesc(log10(image));
title('log10 scale')
colormap gray;
colorbar();
axis tight
hold off

%displaying linear format
figure;
hold on
imagesc((image));
title('linear scale'); 
colormap gray;
colorbar();
axis tight
hold off


shutdown_camera(vid_object);
clear src;
clear vid_object;
