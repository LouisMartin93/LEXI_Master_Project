% Exposure Parameters
num_frames = 100;

video_object = init_camera();
start_camera(video_object);

for i=1:num_frames
    data = take_exposure(video_object);
    phase = process_image_zernike(data);
    phase = phase - phase_ref;
    
    % Use a reconstructor
    control_vector = phase;
    
    %send_new_ref_vec
    
end

stop_camera(video_object);
shutdown_camera(video_object);