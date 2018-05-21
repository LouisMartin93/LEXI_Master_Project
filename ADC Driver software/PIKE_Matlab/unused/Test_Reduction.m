dir = 'Test Images/Right Patch/';


for i = 1:50
    image = take_exposure(vid);
    image = im2double(image);
    fitswrite(image,[dir,'test',num2str(i),'.fits'])
end

