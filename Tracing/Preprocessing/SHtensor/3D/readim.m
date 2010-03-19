% filename = 'C:\Users\arun\Research\TensorVoting\SHtensor\Diadem\DS3_Stack\DS3_02_cropped.tif';
% filename = 'C:\Users\arun\Research\TensorVoting\SHtensor\Diadem\DS5_Stack\OP_3_cropped.tif';
% filename = 'C:\Users\arun\Research\TensorVoting\SHtensor\Diadem\DS5_Stack\OP_2.tif';
function im = readim(filename)
imf = imfinfo(filename);
im1 = imread(filename);
im = zeros([size(im1) size(imf,1) ],'uint8');
for co = 1:size(imf,1)
    im(:,:,co) = imread(filename,co);
end
end