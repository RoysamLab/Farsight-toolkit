im = imread('ecm1.tif');
% im = 0.2989*im(:,:,1)+0.5870*im(:,:,2)+0.1140*im(:,:,3);
imshow(im);
im = double(im);
im1 = im;
min1 = min(im(:));
mean1 = mean(im(:));
max1 =  0.5*max(im(:));
mask = 255*(im>max1);
se = strel('ball',9,9);
mask = imdilate(mask,se);
im = mean1*double(mask/255.0)+double(1-mask/255).*im;
im = uint16(im);
opt.BlackWhite = false;
opt.FrangiScaleRange = [7 9];
opt.FrangiScaleRatio = 1.0;
opt.FrangiBetaOne = 0.5;
opt.FrangiBetaTwo = 15;
[outim,scale,dir] = FrangiFilter2D(double(im),opt);

subplot(2,2,1);imagesc(outim);
subplot(2,2,2);imagesc(dir);
subplot(2,2,3);imagesc(im);
subplot(2,2,4);imagesc(scale);

colorim = zeros([size(outim),3],'uint8');
colorim(:,:,1) = outim*100;
colorim(:,:,2) = im1/3;

figure,imagesc(colorim);
% dir1 = mod(dir+2*pi,pi);
% hist(dir1(:));

% step = 10;

thresh = 0.5;
mask1 = outim>thresh & (255-mask)~=0;
[X,Y] = ind2sub(size(outim),find(mask1));
% [X,Y] = meshgrid(1:step:size(im,1),1:step:size(im,2));
U = scale(mask1).*sin(dir1(mask1));
V = scale(mask1).*cos(dir1(mask1));
hold on;

quiver(Y,X,V(:),U(:));

