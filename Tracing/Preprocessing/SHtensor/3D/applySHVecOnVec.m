function result = applySHTV(TVmodel,gradimg)


L = TVmodel.L;

fftwplanner = 1;



s = 0;

szinput = size(gradimg);
szinput = szinput(2:4);


NormDs = STmultRealF(gradimg,gradimg,0);
NormDs = sqrt(abs(NormDs(1,:,:,:)));
epsi = max(NormDs(:))*0.01;
NormD = split2inter(NormDs,0);
gradimg(1,:,:,:) = gradimg(1,:,:,:) ./ (NormDs+epsi);
gradimg(2,:,:,:) = gradimg(2,:,:,:) ./ (NormDs+epsi);
gradimg(3,:,:,:) = gradimg(3,:,:,:) ./ (NormDs+epsi);
gradimg(4,:,:,:) = gradimg(4,:,:,:) ./ (NormDs+epsi);


display('padding and ffting of basis images');
tic
s1 = TVmodel.bsz(1); s2 = TVmodel.bsz(2); s3 = TVmodel.bsz(3); 
for l = 1:1:L,
    for k = [-1 1],
        TVbasisCUR = myREAL(zeros(2*(l+k+1),szinput(1),szinput(2),szinput(3)));
        p1 = round((szinput(1) - s1)/2);
        p2 = round((szinput(2) - s2)/2);    
        p3 = round((szinput(3) - s3)/2);    
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
result = zeros(4,szinput(1),szinput(2),szinput(3));
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


figure(11);
clf; show(Mag,0.5);


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
    
    







