function writeim(im,filename)
delete(filename);
im  = double(im);
im = im./max(im(:))*255.0;
for co = 1:size(im,3)
    imwrite(uint8(squeeze(im(:,:,co))),filename,'WriteMode','Append','Compression','none');
end
end