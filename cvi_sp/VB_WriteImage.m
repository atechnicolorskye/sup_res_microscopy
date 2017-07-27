function filename = VB_WriteImage(data, folder, tif_name, numFrame)

%% write image
%CCD_imgs =  simparm.y;
filename = [folder, tif_name];     
%numFrame = N; 
% delete previous images
delete(filename);

for nn = 1:numFrame
    imwrite(uint16(data(:,:,nn)),filename,'tif','Compression','none','writemode','append');
end 
