function result = applySHTV(TVmodel)


L = TVmodel.L;

Q = 64;
fftwplanner = 1;


inputimg = (zeros(Q,Q,Q));
[X Y Z] = ndgrid(-(Q/2-1):(Q/2));
%inputimg =   (Z < 5*sin(X/Q*2*pi))+ randn(size(inputimg))*0.1;
inputimg = atan(((X.^2 + Y.^2 + Z.^2) - 15.5^2)*0.001)+ randn(size(inputimg))*0.0;

inputimg = myREAL(inputimg);

s = 0;

szinput = size(inputimg);

epsi = 0.000000000000001;
rs = split2inter(reshape(inputimg,[1 size(inputimg)]),0);
% gradimg = single(get_vesselness_in_spherical(inputimg));
gradimg = STderivUp(rs);
NormDs = STmultRealF(gradimg,gradimg,0);
NormDs = sqrt(abs(NormDs(1,:,:,:)));
NormD = split2inter(NormDs,0);
NormDinv = split2inter(1./(NormDs+epsi),0);
gradimg = STmultRealF(NormDinv,gradimg,1);


NormD = NormD *0 ;

NormD(1,16,32,32) = 1;
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


NormD = circshift(NormD,[0 s s s]);

%figure(1);
%clf; imagesc(real(reshape(NormD(1,:,:,32),[Q Q])));
%drawnow;

%mask = NormD > 0.0 ;% & rand(size(NormD))>0.5;

display('padding and ffting of basis images');
tic
s1 = TVmodel.bsz(1); s2 = TVmodel.bsz(2); s3 = TVmodel.bsz(3); 
for l = 1:1:L,
    ['L=' int2str(l)]
    for k = [-1 1],
        TVbasisCUR = myREAL(zeros(2*(l+k+1),size(inputimg,1),size(inputimg,2),size(inputimg,3)));
        p1 = round((size(inputimg,1) - s1)/2);
        p2 = round((size(inputimg,2) - s2)/2);    
        p3 = round((size(inputimg,3) - s3)/2);    
        TVbasisCUR(:,p1:(p1+s1-1),p2:(p2+s2-1),p3:(p3+s3-1) ) = TVmodel.basis{l,k+2};
        TVbasisCUR = fftshift(fftshift(fftshift(TVbasisCUR,2),3),4);
        TVbasis{l,k+2} = myfftw(TVbasisCUR,[2 3 4],fftwplanner);
    end;
end;
toc;

display('computing evidence images');
tic
E{1} = NormD ;
for l = 2:1:(L+1),   
    E{l} = STmultRealF(E{l-1},gradimg,l-1);      
end;
toc

display('computing fft of e-image');
tic
for l = 1:1:(L+1),         
   Ecur = myfftw(E{l} , [2 3 4],fftwplanner);   
   E{l} = Ecur;
     
end;
toc

display('rendering voting image');
tic
result = zeros(4,size(inputimg,1),size(inputimg,2),size(inputimg,3));
for l = 1:1:L,   
    for k = [-1 1],      
        result = result + STmultRealConvAsymF(TVbasis{l,k+2},E{l+1},1);
    end;
end;

result = myifftw(result,[2 3 4],fftwplanner) / (prod(szinput));
result = circshift(result,[0 1 1 1]);

toc

recon3dZ = reshape(result(end-1,:,:,:),szinput);
recon3dY = -reshape(result(end-2,:,:,:),szinput)*sqrt(2);
recon3dX = -reshape(result(end-3,:,:,:),szinput)*sqrt(2);

figure(9);
Mag = sqrt(recon3dZ.^2 + recon3dX.^2 + recon3dY.^2 );
%clf; slice(double(Mag),q2,q2,q2);
subplot(2,1,1);
clf; orthoview(Mag/max(Mag(:)));




figure(10);
clf;
idx1 = find(X>0 & Y > 0);
[dummy idx] = sort(Mag(:),'descend');
idx = idx(1:500);
quiver3(X(idx),Y(idx),Z(idx),recon3dX(idx),recon3dY(idx),recon3dZ(idx));
axis equal
axis vis3d

figure(11);
clf; show(Mag,8);

%{
figure(1);
clf; slice(double(real(result)),q2,q2,q2);
figure(2);
clf; imagesc(double(reshape(real(result(:,:,q2)),[Q Q])));
figure(3);
clf; show(double(real(result)),1);
%}


return
















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
    
    







