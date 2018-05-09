function Set_Shutter(video_object,Shutter_time)
    source = getselectedsource(video_object);
    source.Shutter = Shutter_time;
end