function steerableTVbasis = expandVF(L)
steerableTVbasis = [];

%%%%%%%%%%%%%%%%% compute basic 2d-voting field

M = 20; %% size of orig vf
figure(1);
sig = 0.3;
%sig = 0.2;
Q = computevec([1 0],M,M,sig,sig*0.6); quiver(real(Q),imag(Q))
%Q = ifft2(fft2(Q) .* fftshift(fspecial('gaussian',size(Q),10)));
Q=Q';
VotefunZ = real(Q);
VotefunW = imag(Q);



%%%%%%%%%%%%%%%%%%%% project 2d-voting field on Legendre basis


%L =7; % degree of expansion
N = 256; %% number of steps on halfcircle
R = 256; %% number of radii

VFsize = 0.08;

%VFsize = 0.08;
%VFsize = 0.03;

expitheta =  exp(pi*i*(1:N)/N);
Radius = ((0:(R-1))'/R*(M-1)/2);
interpolPts = (M+2)/2*(1+i) + Radius * expitheta;
Gauss = exp(-Radius.^2*0.01^2);
polImgZ = interp2(VotefunZ,real(interpolPts),imag(interpolPts),'bilinear');
polImgW = interp2(VotefunW,real(interpolPts),imag(interpolPts),'bilinear');

figure(1);
clf;
quiver(VotefunZ,VotefunW);
%subplot(2,1,1);
%imagesc(VotefunZ);
%subplot(2,1,2);
%imagesc(VotefunW);
%hold on
%plot(real(interpolPts),imag(interpolPts),'*')
figure(2);
subplot(2,1,1);
imagesc(polImgZ)
subplot(2,1,2);
imagesc(polImgW)

costheta = real(expitheta);
sintheta = sqrt(1-costheta.^2);
reconZ = 0;
reconW = 0;

for l = 0:1:(L+1),
    P = legendre(l,costheta)/sqrt(N/pi);
    LegExpansionZ(:,l+1) = polImgZ * (sintheta.*P(1,:))';    
    reconZ = reconZ + LegExpansionZ(:,l+1) * P(1,:) *(2*l+1);
end;
for l = 1:1:(L+1),
    P = legendre(l,costheta)/sqrt(N/pi);
    LegExpansionW(:,l+1) = polImgW * (sintheta.*P(2,:))' ;    
    reconW = reconW + LegExpansionW(:,l+1) * P(2,:)*(2*l+1);
end;

for l = 1:L,
    for k = -1:1,
     % Zproj(:,l,k+2) = LegExpansionZ(:,l+k+1) *ClebschGordan(l,l+k,1,0,0,0);
     % Zproj(:,l,k+2) = Zproj(:,l,k+2)  - 1/sqrt(l*(l+1))* LegExpansionW(:,l+k+1)*ClebschGordan(l,l+k,1,0,1,1);
       Zproj(:,l,k+2) = LegExpansionZ(:,l+k+1)*ClebschGordan(l,l+k,1,0,0,0) - 1/sqrt(l*(l+1))*LegExpansionW(:,l+k+1)*ClebschGordan(l,l+k,1,0,1,1);
      Zproj(:,l,k+2) = Gauss .* Zproj(:,l,k+2) * (2*l+1)*(2*(l+k)+1);
      
    end;
end;
%{      
figure(4);
subplot(2,1,1);
imagesc(LegExpansionZ);
subplot(2,1,2);
imagesc(Zproj(:,:,3));

figure(3);
subplot(2,1,1);
imagesc(reconZ)
subplot(2,1,2);
imagesc(reconW)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%% compute steerable kernels

CS = M+1; %% size of kernel

[X Y Z] = ndgrid(((-(CS-1)/2):((CS-1)/2)));
D = sqrt(X.^2 + Y.^2 + Z.^2)+0.000000000000000001;
gauss = reshape(exp(-D.^2/10),[1 CS CS CS]);
Dinteger = round(D*R*VFsize)+1;
Dvalid = find(Dinteger <= R);
Dmap = Dinteger(Dvalid);

W = (X+i*Y);
SHB(1,:,:,:) = -(W./D)/sqrt(2) ;%.* exp(-D.^2*0.01);
SHB(2,:,:,:) = Z./D ;% .* exp(-D.^2*0.01);
SHB = myREAL(SHB);

SHB = real(split2inter(real(SHB),imag(SHB)));

tic
SH{1} = myREAL(ones(2,CS,CS,CS));
SH{1}(2,:,:,:) = 0;
SH{2} = SHB;
for l = 3:1:(L+2),   
    SH{l} =  STmultRealF(SH{l-1},SHB,l-1) ;   
end;
toc


%clf; showSubband(X,Y,Z,SH{2},0.001),

recon3dZ = 0;
recon3dY = 0;
recon3dX = 0;


tic
for l = 1:1:(L),
    for k = [-1 1 ],
        Rdep = zeros(1,CS,CS,CS);
        Rdep(1,Dvalid) = Zproj(Dmap,l,k+2);
        for m = 1:size( SH{l+k+1},1),
            Basis{l,k+2}(m,:,:,:) = SH{l+k+1}(m,:,:,:) .* Rdep;
        end;
        fprintf('l=%i, k=%i   ->  %f \n',l,k,sum(abs(SH{l+k+1}(:))));
    end;
end;
toc



%%%%%%%% test

% % for l = 1:1:(L+1),     
% %     SHtest{l} = SH{l}(1:2:end,:,:,:) + i* SH{l}(2:2:end,:,:,:) ;
% % end;
% % 
% % for l = 1:1:(L-1),
% %     for k = [-1 0 1 ],
% %         TBasis{l,k+2} = Basis{l,k+2}(1:2:end,:,:,:) + i *Basis{l,k+2}(2:2:end,:,:,:);
% %     end;
% % end;
% % 
% % 
% % recon3d = 0;
% % for l = 1:1:(L-1),
% %     for k = [-1 1 ],
% %        recon3d = recon3d + STmultRealAsymF(TBasis{l,k+2},SHtest{l+1}(:,:,:,:),1);
% %     end;
% % end;
% % 
% % 
% % recon3dZa = real(reshape(recon3d(2,:,:,:),[CS CS CS]));
% % recon3dXa = -real(reshape(recon3d(1,:,:,:),[CS CS CS]))*sqrt(2);
% % recon3dYa = -imag(reshape(recon3d(1,:,:,:),[CS CS CS]))*sqrt(2);
% % 
% % recon3dA(1,:,:,:) = recon3dXa;
% % recon3dA(2,:,:,:) = recon3dYa;
% % recon3dA(3,:,:,:) = recon3dZa;
% % 
% % 
% % cu = 20:30;
% % recon3dX = recon3dXa(cu,cu,cu);
% % recon3dY = recon3dYa(cu,cu,cu);
% % recon3dZ = recon3dZa(cu,cu,cu);
% % 
% % X = X(cu,cu,cu);
% % Y = Y(cu,cu,cu);
% % Z = Z(cu,cu,cu);
% % 
% % figure(9);
% % clf
% % Mag = sqrt(recon3dZ.^2 + recon3dX.^2 + recon3dY.^2 );
% % [dummy idx] = sort(Mag(:),'descend');
% % %idx = idx(1:500);
% % %quiver3(X(idx),Y(idx),Z(idx),recon3dX(idx),recon3dY(idx),recon3dZ(idx));
% % %clf; show(X,Y,Z,Mag,10)
% % %orthoview(Mag/max(Mag(:)));
% % 
% % qq = 250;
% % recon3dX(idx(qq:end)) = 0;
% % recon3dY(idx(qq:end)) = 0;
% % recon3dZ(idx(qq:end)) = 0;
% % tt = 1000;
% % recon3dX = recon3dX ./ (Mag+tt);
% % recon3dY = recon3dY ./ (Mag+tt);
% % recon3dZ = recon3dZ ./ (Mag+tt);
% % 
% % 
% % hc = coneplot(X,Y,Z,-recon3dX,-recon3dY,-recon3dZ,0.1,'nointerp');
% % set(hc,'FaceColor','red','EdgeColor','none');
% % camlight left; lighting phong;
% % view(213,56);
% % 
% % return;
% % 
% % figure(10);
% % clf;
% % Mag = Mag / max(Mag(:));
% % imagesc(reshape(Mag(:,25,:),[49 49]));
% % hold on;
% % X2 = meshgrid(1:2:49);
% % quiver(X2,X2',reshape(recon3dZ(1:2:end,25,1:2:end),[25 25]), reshape(recon3dX(1:2:end,25,1:2:end),[25 25]),1,'y')
% % axis equal
% % axis off;
% % 




%%%%%%%%%%%% test ende

%{
figure(5);
half = (CS+1)/2;
subplot(2,2,1);
imagesc(real(reshape(recon3dZ(:,half,:),[CS CS])));
subplot(2,2,2);
imagesc(real(reshape(recon3dZ(half,:,:),[CS CS])));

figure(6);
clf; show(X,Y,Z,recon3dZ,6),

figure(7);
half = (CS+1)/2;
subplot(2,2,1);
imagesc(real(reshape(recon3dX(half,:,:),[CS CS])));
subplot(2,2,2);
imagesc(real(reshape(recon3dY(half,:,:),[CS CS])));

figure(8);
clf; show(X,Y,Z,Mag,15),
%}

steerableTVbasis.L = L;
steerableTVbasis.basis =  Basis;
steerableTVbasis.bsz = [CS CS CS];



function show(X,Y,Z,gt,thres)

    p = patch(isosurface(X,Y,Z,double(real(gt)),thres));    
    hold on;
    p2 = patch(isosurface(X,Y,Z,-double(real(gt)),thres));
    hold off;

    %isonormals(real(gt),p);
    set(p,'FaceColor','red','EdgeColor','none');
    set(p2,'FaceColor','blue','EdgeColor','none');
    daspect([1 1 1])
    view(3); axis tight
    camlight 
    
    lighting flat;
   % axis([-1 1 -1 1 -1 1]*3);
    axis off

function showSubband(X,Y,Z,ST,thres),
sz = size(ST);    
L = sz(1);
ad = sz(2:4);

for m = 1:L,
   subplot(2,L,m,'align');
   show(X,Y,Z,reshape(ST(m,:),ad),thres);                    
   subplot(2,L,m+L,'align');
   imagesc(real(reshape(ST(m,:,:,32),ad(1:2))));                   
end;
    
    






