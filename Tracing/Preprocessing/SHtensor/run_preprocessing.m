%run_preprocessing
tic
farsight_trunk = 'C:\Users\arun\Research\Farsight\src\trunk';
addpath([ farsight_trunk '\Tracing\Preprocessing\SHtensor\2D']);
addpath([ farsight_trunk '\Tracing\Preprocessing\SHtensor\3D']);
addpath([ farsight_trunk '\Tracing\Preprocessing\CurveletFiltering\fdct_wrapping_matlab']);
addpath('C:\Users\arun\Research\TensorVoting\homomorphic_filtering');
addpath('C:\Users\arun\Research\TensorVoting\focus_detection');

f0 = 'C:\Users\arun\Research\TensorVoting\SHtensor\Diadem\';
f1 = 'hippocampal_part2_cropped';
im = readim([f0 f1 '.tif']);

[focussed, depthmap] = focus_detect(im);
writeim(focussed,['new_best_focus_' f1 '.tif']);
writeim_norescale(depthmap,['new_depthmap_' f1 '.tif']);
return
focussed = imread(['best_focus_' f1 '.tif']);
h = fspecial('gaussian',150,60);
im2 = filter2(h,(focussed));
writeim(im2,'gaussblur.tif');
focussed = single((focussed - uint8(im2)).*uint8(focussed > im2));
% h1 = fspecial('disk',250);
% varim = sqrt(conv2fft(focussed.^2,h1,'same')-conv2fft(focussed,h1,'same').^2);
% varim = varim./(max(varim(:)))/10;
% figure, imagesc(varim);
% focussed = (focussed > 255*graythresh(focussed)*0.8);
% focussed = uint8(255*focussed);
figure,imshow(uint8(focussed))
[sx,sy] = ginput(1);
writeim(focussed,['corrected_focussed_' f1 '.tif']);

% pause
padding_width = 40;
focussed = add_padding(focussed, padding_width);
% varim = add_padding(varim,padding_width);
imshow(uint8(focussed));
drawnow
[out1,gradmag1] = scalar_voting_main((focussed),5,0.35*255,8);
[out2,gradmag2] = scalar_voting_main((focussed),5,0.15*255,8);
gauss_sig=500;
[x1,y1] = meshgrid(1:size(gradmag1,2),1:size(gradmag1,1));
% immaskh = double(((x1-sx).^2+(y1-sy).^2)<=gauss_sig^2);
immaskh = exp(-((x1-sx).^2+(y1-sy).^2)/2/gauss_sig^2);
immaskl = 1-immaskh;
figure,imshow(immaskh);
figure,imshow(immaskl);
figure,imagesc(immaskh.*double(gradmag1)+immaskl.*double(gradmag2));
gradmag = immaskh.*double(gradmag1)+immaskl.*double(gradmag2);
out = immaskh.*out1+immaskl.*out2;
% figure, imagesc(gradmag1);
% figure, imagesc(gradmag2);

% [out,gradmag] = scalar_voting_main(double(focussed));
gradmag = double(gradmag);
opt.BlackWhite = false;
opt.FrangiScaleRange = [1 5];
% out = out*255./max(out(:));

% out = add_padding(remove_padding(out,3),padding_width+3);
% gradmag = add_padding(remove_padding(gradmag,3),padding_width+3);
% focussed = double(add_padding(focussed,padding_width));

[vesselness_out, scale, dir] = FrangiFilter2D(out*255./max(out(:)),opt);
[vesselness_grad, scale, dir] = FrangiFilter2D(double(gradmag),opt);
% [vesselness4, scale, dir] = FrangiFilter2D(double(gradmag1),opt);
% [vesselness3, scale, dir] = FrangiFilter2D(double(focussed), opt);


out = remove_padding(out,padding_width);
gradmag = remove_padding(gradmag,padding_width);
focussed = remove_padding(focussed,padding_width);
vesselness_out = remove_padding(vesselness_out,padding_width);
vesselness_grad = remove_padding(vesselness_grad,padding_width);


% figure, imagesc(vesselness);
figure, imagesc(uint8(gradmag));
figure, imagesc(vesselness_out);
figure, imagesc(vesselness_grad);

writeim(vesselness_out, ['vesselness_scalar_voting_' f1 '.tif']);
writeim(vesselness_grad,['vesselness_curvelets_' f1 '.tif']);
% writeim(vesselness3,['vesselness_best_focussed_' f1 '.tif']);

writeim(gradmag,['curvelet_output_' f1 '.tif']);
writeim(out,['scalar_voting_output_' f1 '.tif']);

toc
1;