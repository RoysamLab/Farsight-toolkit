function im3D(Vol, fname)
figure(33); 
d = size(Vol);
mx = max(Vol(:));
mn = min(Vol(:));
for i = 1:d(3)
    figure(33); 
    imagesc(Vol(:,:,i),[mn mx]); 
    colormap gray;
%     colormap jet;
    colorbar; axis equal;
    title([num2str(i),'/',num2str(d(3))]);
    if nargin==2
        F(i) = getframe;
    end
    pause(0.5);
end
if nargin==2
    movie2avi(F,fname);
end
    