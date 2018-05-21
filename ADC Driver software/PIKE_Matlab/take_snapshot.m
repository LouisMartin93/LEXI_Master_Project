function image = take_snapshot(video_object)
    
    if isrunning(video_object)
        disp('taking image');
        image = getsnapshot(video_object);
    end
end
