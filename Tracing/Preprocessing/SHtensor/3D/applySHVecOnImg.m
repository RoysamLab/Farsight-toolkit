function [result,Mag] = applySHVecOnImg(TVmodel,img)


L = TVmodel.L;

fftwplanner = 2;


inputimg = myREAL(img);

s = 0;

szinput = size(inputimg);

epsi = 0.000000000000001;
rs = split2inter(reshape(inputimg,[1 size(inputimg)]),0);
% gradimg = STderivUp(rs);
[gradimg,J] = get_vesselness_in_spherical(inputimg);
% [gradimg,J] = get_weingarten_in_spherical(inputimg);
gradimg = single(gradimg);
NormDs = STmultRealF(gradimg,gradimg,0);
NormDs = sqrt(abs(NormDs(1,:,:,:)));
NormD = split2inter(NormDs,0);
NormDinv = split2inter(1./(NormDs+epsi),0);
gradimg = STmultRealF(NormDinv,gradimg,1);


display('computing evidence images');
tic
E{1} = NormD ;
for l = 2:1:(L+1),   
    E{l} = STmultRealF(E{l-1},gradimg,l-1);      
end;
toc

display('computing fft of e-image');
tic
for l = 2:2:(L+1),         
    ['l=' int2str(l)]
   Ecur = myfftw(E{l} , [2 3 4],fftwplanner);   
   E{l} = Ecur;
     
end;
toc

% display('padding and ffting of basis images');
% tic
%
% for l = 1:2:L,
%     for k = [-1 1],
% 
%     end;
% end;
% toc;


display('computing ffts of basis and rendering voting image');
tic
s1 = TVmodel.bsz(1); s2 = TVmodel.bsz(2); s3 = TVmodel.bsz(3); 
result = zeros(4,size(inputimg,1),size(inputimg,2),size(inputimg,3));
for l = 1:2:L,   
    for k = [-1 1],      
        ['L = ' int2str(l) ' k = ' int2str(k)]
        TVbasisCUR = myREAL(zeros(2*(l+k+1),size(inputimg,1),size(inputimg,2),size(inputimg,3)));
        p1 = round((size(inputimg,1) - s1)/2);
        p2 = round((size(inputimg,2) - s2)/2);    
        p3 = round((size(inputimg,3) - s3)/2);    
        TVbasisCUR(:,p1:(p1+s1-1),p2:(p2+s2-1),p3:(p3+s3-1) ) = TVmodel.basis{l,k+2};
        TVbasisCUR = fftshift(fftshift(fftshift(TVbasisCUR,2),3),4);
        TVbasisCUR = myfftw(TVbasisCUR,[2 3 4],fftwplanner);
        result = result + STmultRealConvAsymF(TVbasisCUR,E{l+1},1);
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
figure(10);clf;
h = vol3d('cdata',Mag,'texture','3D');
vol3d(h);
figure(11);clf;
h1 = vol3d('cdata',img,'texture','3D');
vol3d(h1);
figure(12);clf;
h2 = vol3d('cdata',J,'texture','3D');
vol3d(h2);
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
    
    







