function steerableTVbasis = expandVF

%%%%%%%%%%%%%%%%% compute basic 2d-voting field

M = 30; %% size of orig vf
figure(1);
sig = 0.15;
sig2 = sig*3.0;
Q = computevec([1 0],M,M,sig,sig2);%quiver(real(Q),imag(Q));
  Q = Q';
figure(50);imagesc(abs(Q));

Votefun = abs(Q);
hold on;
quiver(real(Q),imag(Q));
% figure
%%%%%%%%%%%%%%%%%%%% project 2d-voting field on Legendre basis


L = 8; % degree of expansion
N = 256; %% number of steps on halfcircle
R = 256; %% number of radii

VFsize = 0.08;

expitheta =  exp(pi*i*(1:N)/N);
interpolPts = (M+2)/2*(1+i) + ((0:(R-1))'/R*(M-1)/2) * expitheta;

polImg = interp2(Votefun,real(interpolPts),imag(interpolPts),'bilinear');
figure(1),imagesc(Votefun);
hold on
plot(real(interpolPts),imag(interpolPts),'*')
figure(2);
imagesc(polImg)

costheta = real(expitheta);
sintheta = sqrt(1-costheta.^2);
LegExpansion = zeros(R,L);
recon = 0;

for l = 0:2:L,
    P = legendre(l,costheta,'sch')/sqrt(N/pi);
    LegExpansion(:,l/2+1) = polImg * (sintheta.*P(1,:))' * (2*l+1);  
    recon = recon + LegExpansion(:,l/2+1) * P(1,:);
end;
figure(4);
imagesc(LegExpansion);
figure(3);
imagesc(recon)



%%%%%%%%%%%%%%%%%%%%%%%%%%% compute steerable kernels

CS = M+1; %% size of kernel

[X Y Z] = ndgrid(((-(CS-1)/2):((CS-1)/2)));
D = sqrt(X.^2 + Y.^2 + Z.^2)+0.00000001;
Dinteger = round(D*R*VFsize)+1;
Dvalid = find(Dinteger <= R);
Dmap = Dinteger(Dvalid);

W = (X+i*Y);
SHB(1,:,:,:) = -(W./D)/sqrt(2) ;%.* exp(-D.^2*0.01);
SHB(2,:,:,:) = Z./D ;% .* exp(-D.^2*0.01);
SHB = myREAL(SHB);

SHB = real(split2inter(real(SHB),imag(SHB)));

% clf; showSubband(X,Y,Z,SHB,0.001),

tic
SH{1} = myREAL(ones(2,CS,CS,CS));
SH{1}(2,:,:,:) = 0;
SH{2} = STmultiply(SHB,SHB,2);
for l = 4:2:L,   
    SH{l/2+1} =  STmultiply(SH{2},SH{l/2},l) ;
end;
toc


% clf; showSubband(X,Y,Z,SH{2},0.001),

recon3d = 0;
tic
for l = 0:2:L,
    Rdep = zeros(1,CS,CS,CS);
    Rdep(1,Dvalid) = LegExpansion(Dmap,l/2+1);
    for m = 1:size(SH{l/2+1},1);   
        SH{l/2+1}(m,:,:,:) = SH{l/2+1}(m,:,:,:) .* Rdep;
    end;
    recon3d = recon3d + reshape(SH{l/2+1}(end-1,:,:,:),[CS CS CS]) ;
end;
toc

figure(5);
half = (CS+1)/2;
subplot(2,2,1);
imagesc(real(reshape(recon3d(:,half,:),[CS CS])));
subplot(2,2,2);
imagesc(real(reshape(recon3d(half,:,:),[CS CS])));
subplot(2,2,3);
imagesc(Votefun);
subplot(2,2,4);
imagesc(real(reshape(recon3d(:,:,half),[CS CS])));
figure(6);
clf; show(X,Y,Z,recon3d,2),

figure(7);
subplot(2,2,1);
imagesc(squeeze(max(real(recon3d),[],1)));
subplot(2,2,2);
imagesc(squeeze(max(real(recon3d),[],2)));
subplot(2,2,3);
imagesc(squeeze(max(real(recon3d),[],3)));

steerableTVbasis.L = L;
steerableTVbasis.basis = SH;
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

function showSubband(X,Y,Z,ST,thres)
sz = size(ST);    
L = sz(1);
ad = sz(2:4);

for m = 1:L,
   subplot(2,L,m,'align');
   show(X,Y,Z,reshape(ST(m,:),ad),thres);                    
   subplot(2,L,m+L,'align');
   imagesc(real(reshape(ST(m,:,:,25),ad(1:2))));                   
end;
    
    






