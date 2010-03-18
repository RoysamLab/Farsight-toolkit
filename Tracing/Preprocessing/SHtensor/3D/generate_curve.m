

t = -10:0.01:10;

x = t;
y = sin(t/8);
z = t.^5;
% 
x = 5 + (t+10).*cos(t*pi/5);
y = 5 + (t+10).*sin(t*pi/5);
z = t;
% xtan = xtan ./ normtan;
% ytan = ytan ./ normtan;
% ztan = ztan ./ normtan;


% 
% plot3(x,y,z);
xsize = 200;
ysize = 200;
zsize = 50;
im = zeros(xsize,ysize,zsize,'uint8');

numt = numel(t);
xmin = min(x);
ymin = min(y);
zmin = min(z);
xmax = max(x);
ymax = max(y);
zmax = max(z);
w = 1;
x1 = x;
y1 = y;
z1 = z;
for co = 1:numel(t)
    if( mod(co,100) <20)
        continue;
    end
%     if( co <100 ||  co > 400)
%         continue;
%     end
    x1(co) = int32((x(co)-xmin)/(xmax-xmin)*(xsize-1)+ 1);
    y1(co) = int32((y(co)-ymin)/(ymax-ymin)*(ysize-1)+ 1);
    z1(co) = int32((z(co)-zmin)/(zmax-zmin)*(zsize-1)+ 1);
    x2 = [max(1,x1(co)-w):min(xsize,x1(co)+w)];
    y2 = [max(1,y1(co)-w):min(ysize,y1(co)+w)];
    z2 = [max(1,z1(co)-w):min(zsize,z1(co)+w)];
    im( x2,y2,z2) = 127;
    x1(co) = x1(co) + 10;
    y1(co) = y1(co) + 10;
    z1(co) = z1(co) + 10;
    x2 = [max(1,x1(co)-w):min(xsize,x1(co)+w)];
    y2 = [max(1,y1(co)-w):min(ysize,y1(co)+w)];
    z2 = [max(1,z1(co)-w):min(zsize,z1(co)+w)];
%     im( x2,y2,z2) = 255;
end
width = 3;
coords = [50,50,50;
           9,22,31;
            93,82,74];
% for co = 1:size(coords,1)
%     im(coords(co,1)-width:coords(co,1)+width,coords(co,2)-width:coords(co,2)+width,coords(co,3)-width:coords(co,3)+width)=0;
% end
rand1 = rand(size(im))*150;
% im = im + uint8(rand1);
% imshow(squeeze(max(im,[],2)));
h = vol3d('cdata',im,'texture','3d');
vol3d(h);
% opt.BlackWhite = false;
% opt.FrangiScaleRange = [1 7];
% [J,Scale,Vx,Vy,Vz] = FrangiFilter3D(im,opt);

% for co = 1:numel(t)
%     xtan(co) = Vx(int32(x1(co)),int32(y1(co)),int32(z1(co)));
%     ytan(co) = Vy(int32(x1(co)),int32(y1(co)),int32(z1(co)));
%     ztan(co) = Vz(int32(x1(co)),int32(y1(co)),int32(z1(co)));
% end

% figure, h3 = vol3d('cdata',J,'texture','3D');
% vol3d(h3);
%  step1 = 10;
%  figure(1);
%  plot3(x1(1:step1:end),y1(1:step1:end),z1(1:step1:end),'*');
%  hold on;
%  quiver3(x1(1:step1:end),y1(1:step1:end),z1(1:step1:end),xtan(1:step1:end),ytan(1:step1:end),ztan(1:step1:end));
%  figure(2);
%  imagesc(squeeze(max(J,[],3)))
%  figure(3);
%  imagesc(squeeze(max(Scale,[],3)))