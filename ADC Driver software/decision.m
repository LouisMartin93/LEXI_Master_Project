function decision(videoObject)

choice = input('Is the image aligned enough (y or n): ','s');
if choice == 'n'
    cont_exposures(50,videoObject,'log');
    close;
    decision(videoObject)
end

end