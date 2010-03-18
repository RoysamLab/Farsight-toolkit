
im = imread('sample_from_paper.tif');
 im = double(im);
% min1 = min(im(:));
% mean1 = mean(im(:));
% max1 =  0.5*max(im(:));
% mask = 255*(im>max1);
% se = strel('ball',9,9);
% mask = imdilate(mask,se);
% im = mean1*double(mask/255.0)+double(1-mask/255).*im;
% im = uint16(im*5);

%  im = 0.2989*im(:,:,1)+0.5870*im(:,:,2)+0.1140*im(:,:,3);
% imshow(im);
% im = double(im);
% min1 = min(im(:));
% max1 = 0.2* max(im(:));
% im = min1*double(im>max1)+double(im<max1).*im;
% im = uint16(im);
im1 = im;
%% smooth and find gradient
% im = medfilt2(im,[5 5]);
% hyfilt = fspecial('sobel');
% hxfilt = hyfilt';
% Iy = imfilter(double(im),hyfilt,'replicate');
% Ix = imfilter(double(im),hxfilt,'replicate');


gauss_sigma = 1.0;
maskSize = 4*gauss_sigma;

    gmap = double(fspecial('gaussian', maskSize, gauss_sigma));
   
    gLeft  = [gmap, zeros(maskSize,1)];
    gRight = [zeros(maskSize,1), gmap];
    gDiff = gLeft-gRight;
    DoG.Gx = gDiff(:,1:maskSize);

    gTop    = [gmap; zeros(1,maskSize)];
    gBottom = [zeros(1,maskSize); gmap];
    gDiff = gTop-gBottom;
    DoG.Gy = gDiff(1:maskSize,:);
    
Ix = conv2fft(double(im),DoG.Gx,'same');
Iy = conv2fft(double(im),DoG.Gy,'same');

gradmag = sqrt(Ix.^2+Iy.^2);

%%
opt.FrangiScaleRange = [1 3];
opt.BlackWhite = false;
[outim,scale,dir] = FrangiFilter2D(double(im),opt);
imagesc(outim);
1;
 gradmad = outim;
Ix = outim.*cos(dir);
Iy = outim.*sin(dir);
[X,Y] = meshgrid(1:size(im,2),1:size(im,1));
% imshow(outim);hold on;
% quiver(X,Y,Ix,Iy);
% Ix = Ix ./ (gradmag + 1);
% Iy = Iy ./ (gradmag + 1);
% imshow(uint8(gradmag));

% step = 5;
% [x,y] = meshgrid(1:step:size(im,2),1:step:size(im,1));
% hold on;
% Ixs = Ix(1:step:size(im,1),1:step:size(im,2));
% Iys = Iy(1:step:size(im,1),1:step:size(im,2));
% quiver(x,y,Ixs,Iys);

%% compute the initial tensor

Axx = Ix.*Ix ;%+ 0.001.*Iy.*Iy; % Ixp*Ixp
Axy = Ix.*Iy ;%- 0.001.*Ix.*Iy;
Ayy = Iy.*Iy ;%+ 0.001.*Ix.*Ix;

orientation = mod(0.5*angle(Axx - Ayy + 2*1i*Axy)+pi,pi);
stickness = sqrt((Axx+Ayy).^2 - 4*((Axx.*Ayy)-Axy.^2));
ballness = 0.5*((Axx+Ayy).^2 - stickness);
% figure, imshow(uint8(128*stickness));
%  figure, imshow(uint8(128*ballness));
% figure, imshow(uint8(255*orientation/pi),jet(180))

A0 = Axx + Ayy ;
A2 = Axx - 2*1i*Axy - Ayy;
A2n = Axx + 2*1i*Axy - Ayy;

sigma = 1;
K = double(sigma*3);
[x,y] = meshgrid(-K:1:K,-K:1:K);
norm = sqrt(x.^2+y.^2);
norm(K+1,K+1) = 1;
mat = (x+1i*y)./norm;
mat(K+1,K+1)=1;
for co = 0:2:8
    w{co/2+1} = exp(-(x.^2+y.^2)/2/sigma/sigma).*(mat)^co;
end

for co = 0:2:6
    c{co/2+1} = stickness.*exp(-1i*co*orientation);
end

U2n = conv2fft(conj(c{2}),w{1},'same') + 4.0*conv2fft(c{1},w{2},'same')+6.0*conv2fft(c{2},w{3},'same') + 4*conv2fft(c{3},w{4},'same') + conv2fft(c{4},w{5},'same');
U2 = conj(U2n);
U0 = real(6.0*conv2fft(c{1},w{1},'same') + 8.0*conv2fft(c{2},w{2},'same') + 2*conv2fft(c{3},w{3},'same'));

sticknew = abs(U2n);
orientnew = mod(0.5*angle(U2n)+pi,pi);
ballnew = 0.5*(U0 - abs(U2));
figure
subplot(2,2,1); imagesc(im1);
subplot(2,2,2); imagesc(uint8(gradmag));
subplot(2,2,3); imagesc(uint8(255*sticknew/max(sticknew(:))))
% subplot(2,2,4); imagesc(ballnew);