for i = 1:10
    snapshot = getsnapshot(vidobj);
    imagesc(snapshot);
    colorbar();
end