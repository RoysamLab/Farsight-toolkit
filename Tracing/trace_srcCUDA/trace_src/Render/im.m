function im(img,figno)
if nargin ==1
    figno = 1;
end
figure(figno);
imagesc(img); colormap gray;
colorbar;
axis off normal;

