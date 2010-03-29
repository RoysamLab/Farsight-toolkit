function [result] = applySHTV(TVmodel,im,cness)

% pad it
if mod(size(im,1),2) ~= 0
    pad1 = TVmodel.bsz(1);
else
    pad1 = TVmodel.bsz(1)-1;
end
if mod(size(im,2),2) ~= 0
    pad2 = TVmodel.bsz(2);
else
    pad2 = TVmodel.bsz(2)-1;
end
if mod(size(im,3),2) ~= 0
    pad3 = TVmodel.bsz(3);
else
    pad3 = TVmodel.bsz(3)-1;
end

imnew =  zeros([size(im,1) + pad1, size(im,2) + pad2, size(im,3)+pad3],'single');
imnew(1:size(im,1),1:size(im,2),1:size(im,3)) = im;
im = imnew;
clear imnew
cnessnew = zeros(size(cness,1),size(cness,2)+pad1,size(cness,3)+pad2,size(cness,4)+pad3,'single');
cnessnew(:,1:size(cness,2),1:size(cness,3),1:size(cness,4))=cness; 
cness = cnessnew;
clear cnessnew

L = TVmodel.L;

Q = 64;%size(im,1);
fftwplanner = 1;


inputimg = (zeros(Q,Q,Q));
[X Y Z] = ndgrid(-(Q/2-1):(Q/2));
%inputimg =   (Z < 5*sin(X/Q*2*pi))+ randn(size(inputimg))*0.1;
% inputimg = atan(((X.^2 + Y.^2 + Z.^2) - 5^2)*0.01)+ randn(size(inputimg))*0;
% inputimg(Q/2:end,Q/2:end,Q/2:end)=0;
inputimg = im;
Q = size(im,1);
inputimg = myREAL(inputimg);


s = 0;

szinput = size(inputimg);

% figure(4);
q2 = round(Q/2);
% clf; slice(double(real(inputimg)),q2,q2,q2);

epsi = 0.000000001;
% rs = split2inter(reshape(inputimg,[1 size(inputimg)]),0);
% gradimg = single(get_curvelet_in_spherical(inputimg));
gradimg = single(cness);
% gradimg = STderivUp(rs);
% NormDs = 2*(gradimg(1,:,:,:).^2+gradimg(2,:,:,:).^2)+gradimg(3,:,:,:).^2;
NormDs = STmultiply(gradimg,gradimg,0);
NormDs = sqrt(abs(NormDs(1,:,:,:)));
NormD = split2inter(NormDs,0);
NormDinv = split2inter(1./(NormDs+epsi),0);
gradimg = STmultiply(NormDinv,gradimg,1);
clear NormDs NormDinv cness im
% figure(10);slice(double(real(squeeze(NormD(1,:,:,:)))),q2,q2,q2);
% figure(11);slice(double(real(squeeze(gradimg(1,:,:,:)))),q2,q2,q2);
% figure(12);slice(double(real(squeeze(gradimg(2,:,:,:)))),q2,q2,q2);
% figure(13);slice(double(real(squeeze(gradimg(3,:,:,:)))),q2,q2,q2);
% figure(14);slice(double(real(squeeze(gradimg(4,:,:,:)))),q2,q2,q2);
gradimg;
% NormD(1,Q/2:end,Q/2:end,Q/2:end)=0;
%  NormD = NormD *0 ;
%  NormD(1,50,50,50) =1;
%  NormD(1,9,22,31) = 1;
%  NormD(1,96,89,83) = 1;
%  NormD(1,92,79,70) = 1;
%  NormD(1,26,44,49) = 1;
% 
% NormD(1,16,32,32) = 1;
% NormD(1,32,16,32) = 1;
% NormD(1,48,32,32) = 1;
% NormD(1,32,48,32) = 1;
% NormD(1,32,32,48) = 1;
% NormD(1,32,32,16) = 1;
% 
% NormD(1,21,43,32) = 1;
% NormD(1,43,21,32) = 1;
% NormD(1,21,21,32) = 1;
% NormD(1,43,43,32) = 1;
% 
% NormD(1,21,32,21) = 1;
% NormD(1,43,32,43) = 1;
% NormD(1,21,32,43) = 1;
% NormD(1,43,32,21) = 1;
% 
% NormD(1,32,21,21) = 1;
% NormD(1,32,43,43) = 1;
% NormD(1,32,21,43) = 1;
% NormD(1,32,43,21) = 1;
% 
% 
% NormD = circshift(NormD,[0 s s s]);

%figure(1);
%clf; imagesc(real(reshape(NormD(1,:,:,32),[Q Q])));
%drawnow;

%mask = NormD > 0.0 ;% & rand(size(NormD))>0.5;

display('padding and ffting of basis images');
tic
s1 = TVmodel.bsz(1); s2 = TVmodel.bsz(2); s3 = TVmodel.bsz(3); 
for l = 0:2:L,
    TVbasisCUR = zeros(2*(l+1),size(inputimg,1),size(inputimg,2),size(inputimg,3),'single');
    p1 = round((size(inputimg,1) - s1)/2)+1;
    p2 = round((size(inputimg,2) - s2)/2)+1;    
    p3 = round((size(inputimg,3) - s3)/2)+1;    
    TVbasisCUR(:,p1:(p1+s1-1),p2:(p2+s2-1),p3:(p3+s3-1) ) = TVmodel.basis{l/2+1};
    TVbasisCUR = fftshift(fftshift(fftshift(TVbasisCUR,2),3),4);
    TVbasis{l/2+1} = myfftw(TVbasisCUR,[2 3 4],fftwplanner);
    clear TVbasisCUR;
end;
toc;

display('computing evidence images');
tic
E{1} = NormD ;
E2 = STmultiply(gradimg,gradimg,2);
for l = 2:2:L,   
    E{l/2+1} = STmultiply(E2,E{l/2},l);      
end;
toc

display('computing fft of e-image');
tic
for l = 0:2:L,         
   Ecur = myfftw(E{l/2+1} , [2 3 4],fftwplanner);   
   E{l/2+1} = Ecur;
     
end;
toc

display('rendering voting image');
tic
result = zeros(2,size(inputimg,1),size(inputimg,2),size(inputimg,3));
for l = 0:2:L,   
    result = result + STmultiplyFourier(TVbasis{l/2+1},E{l/2+1},0);
end;

result = myifftw(result,[2 3 4],fftwplanner) / (prod(szinput));

result = reshape(result(1,:,:,:),szinput);
toc

result = result(1:size(result,1)-pad1,1:size(result,2)-pad2,1:size(result,3)-pad3);
%result = smooth3(result);
minr = min(result(:));
maxr = max(result(:));
% result = (result > 2) .*2 +(result<=2).*result;
% result = 10.0*(result - minr)./(maxr-minr);
% figure(1);
% clf; slice(double(real(result)),q2-5,q2+5,q2);
% figure(2);
% clf; imagesc(double(reshape(real(result(:,q2+15,:)),[Q Q])));
% figure(3);
% clf; show(double(real(result)),1);
% imagesc(squeeze(max(result,[],3)));
 result1 = (result>20).*result;
figure, h = vol3d('cdata',result1,'texture','3d');
vol3d(h);
% writeim(result,'result.tif');
% writeim(im,'input.tif');
% writeim(squeeze(NormD(1,:,:,:)),'normD.tif');
% figure
1;
end



















function show(gt,thres)

    p = patch(isosurface(real(gt),thres));
    hold on;
 %   p2 = patch(isosurface(-real(gt),thres));
 %   hold off;

    %isonormals(X,Y,Z,G,p);
    set(p,'FaceColor','red','EdgeColor','none');
  %  set(p2,'FaceColor','blue','EdgeColor','none');
    daspect([1 1 1])
    view(3); axis tight
    camlight 
    lighting gouraud
   % axis([-1 1 -1 1 -1 1]*32 +32);%
    axis([1 size(gt,1) 1 size(gt,2) 1 size(gt,3) ]);
    %axis off

end


function showSubband(ST,thres),
sz = size(ST);    
L = sz(1);
ad = sz(2:4);

for m = 1:L,
   subplot(2,L,m,'align');
   show(reshape(ST(m,:),ad),thres);                    
   subplot(2,L,m+L,'align');
   imagesc(real(reshape(ST(m,:,:,32),ad(1:2))));                   
end;
    
    
end






