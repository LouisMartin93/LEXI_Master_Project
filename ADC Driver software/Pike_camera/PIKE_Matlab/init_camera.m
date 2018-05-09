function [video_object] = init_camera()
%INIT_CAMERA Summary of this function goes here
%   Detailed explanation goes here
    
    % Make a list of connected cameras
    current_hw = imaqhwinfo;
    disp(current_hw.InstalledAdaptors);
    
    % Choose a camera
    camera_type = input('Choose a camera: ', 's');
    
    % Search for Firewire cameras
    dev_info = imaqhwinfo(camera_type, 1);
    
    % Display the video formats
    celldisp(dev_info.SupportedFormats);
    
    frame_format = input('Choose a frame format: ', 's');
    
    % Init camera in a certain video format
    video_object = videoinput(camera_type, 1, frame_format);
    
    % Gets properties of the camera;
    src = getselectedsource(video_object);
    % choose a ROI
    ROI = input('Select camera ROI ([offsetx offsety widthx widthy]):');
    if isempty(ROI)
        ROI = [0 0 640 480]
    end
    video_object.ROIPosition = ROI;
    
    % choose a shutter time
    Shutter = input('Select Shutter time (ms): ');
    src.Shutter = Shutter;
    
    % Choose number of frames per trigger (default is 10)
    trig  = input('Enter number of frames per trigger: ');
    video_object.FramesPerTrigger = trig;
end

