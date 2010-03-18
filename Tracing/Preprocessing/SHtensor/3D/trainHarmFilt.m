function model = train(images,labels,param,regu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  model = train(images,labels,param,regu)
%
%  images - cell array of source images
%  labels - cell array of corresponding label/target images
%  param  - array of parameters [L gin gout gamma]
%             L      - degree of expansion of harmonic function
%             gin    - width of gaussian on source image
%             gout   - width of gaussian to generate target image
%             gamma  - regularization param. of stddev-normalization 
%                      (gamma < 0 => no regularization)                     
%  regu   - regulariztion parameter for solving the lsq-regression
%  model  - trained model for passing to applyHarmFilt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = param(1);
gwidin =param(2);
gwidout = param(3);
gamma = param(4);


if nargin == 3,
    regu = 0;
end;



bo = 0;             %% boarder, just voxels inside the boarder are considered as training examples
fftwplanner = 2;    %% fftw planning flag


FIcol = []; %% features
licol = []; %% targets
nVcol = 0;  %% number of voxels per image

for nI = 1:length(images),
    imgorig = images{nI};       
       
   
    img = myREAL(imgorig);    
    imgsz = size(img);
    
    labelimg = myREAL(labels{nI});
    if isempty(labelimg),
        labelimg = zeros(imgsz);
    end;
    
   
    gaussin = GaussCache(imgsz,gwidin,fftwplanner);
    gaussout = GaussCache(imgsz,gwidout,fftwplanner);
 
               
    display('blurring input');
    tic
    inputimg = myREAL(img);
    inputimg = reshape(inputimg,[1 imgsz]);
    inputimg(2,:,:,:) = 0;      
    inputimg_sav = inputimg;
    inputimg = mysmooth(inputimg,gaussin,fftwplanner)/prod(imgsz);
    toc

    display('computing stddev');
    tic
    simg = mysmooth(inputimg_sav.^2,gaussin,fftwplanner)/prod(imgsz);
    stddev = simg(1,:,:,:) - inputimg(1,:,:,:).^2;
    stddev = real(sqrt(stddev.* (stddev>0)));
    toc
    
    display('computing derivatives');
    tic
    D{1} = inputimg;
    inrank(1) = 0;
    for k = 2:L,
        D{k} = STderivUp(D{k-1});
        if k == 2 && gamma >= 0,         
            normD2 = gamma+stddev;
            D{2}(1,:,:,:) =  D{2}(1,:,:,:) ./normD2;
            D{2}(2,:,:,:) =  D{2}(2,:,:,:) ./normD2;
            D{2}(3,:,:,:) =  D{2}(3,:,:,:) ./normD2;
        end;
        inrank(k) = k-1;
        fprintf('.');
    end;
   % D{L+1}(1,:,:,:) = ones(imgsz);
    %D{L+1}(2,:,:,:) = 0;
    %inrank(L+1) = 0;
  %  D{1}(1,:,:,:) = 1;    
    toc

    %% compute product/feature types [k j l], eq. to D_k *_l D_j
    %% k,j - rank+1 of spherical derivatives
    %% l - rank of result of spherical product
    cnt = 1;
    for k = 1:length(inrank),
        for j = k:length(inrank),
            for l = abs(inrank(k)-inrank(j)):(inrank(k)+inrank(j)),
                if mod(l+inrank(k)+inrank(j),2) == 0 & (l <= L),% &  inrank(k) > 0 & inrank(j) > 0,% &  mod(inrank(k)+inrank(j),2) == 0
                    ranktup(cnt,:) = [k j l];
                    cnt = cnt + 1;
                end;
            end;
        end;
    end;
  
    
    display('computing featimgs');
    tic
    clear FeatImg;
    FeatImg = zeros([size(ranktup,1),imgsz]);
    for k = 1:size(ranktup,1),
        tic
        Fimg = STmultiply(D{ranktup(k,1)},D{ranktup(k,2)},ranktup(k,3));
        for j = 1:ranktup(k,3),
            Fimg = STderivDown(Fimg);
        end;
       
        fprintf('%i/%i multi/deriv  %f\n',k,size(ranktup,1),toc);
        tic
        Fimg = mysmooth(Fimg,gaussout,fftwplanner) / prod(imgsz);
        FeatImg(k,:,:,:) = Fimg(1,:,:,:);
        fprintf('      blur %f \n ',toc);
        
    end;
    toc
    FI = FeatImg(:,(bo+1):(end-bo),(bo+1):(end-bo),(bo+1):(end-bo));
    li = labelimg((bo+1):(end-bo),(bo+1):(end-bo),(bo+1):(end-bo));
       
    FI = double(FI(:,:));
    numVoxels = size(FI,2);

    %% collect training examples 
    FIcol = [FIcol FI];
    licol = [licol ; double(li(:))];
    nVcol = nVcol + numVoxels;
    
end

% normalize features
renormfac = sqrt(sum(FIcol.^2,2)/nVcol);
invrenormfac = 1./renormfac;
FIcol = FIcol.* (invrenormfac*ones(1,nVcol));

% compute covariances
C = FIcol*FIcol';
b = FIcol*double(licol(:));


    
display('solve regression problem');
tic
alpha = (C+regu*diag(invrenormfac.^2))\b;
alpha = alpha / (alpha'*b) *sum(licol);
alpha = invrenormfac.*alpha;
alpha = myREAL(alpha);
toc


% compute result on last training image
result = reshape(alpha'*FeatImg(:,:),imgsz);


% show results
img2show = reshape(inputimg(1,:,:,:),size(result));
figure(2);
result = (result>0).*result;
resimg(:,:,:,1) = result/max(result(:));
resimg(:,:,:,2) = img2show.*(img2show>0) / max(img2show(:));
resimg(:,:,:,3) = 0*imgorig.*(imgorig>0) / max(imgorig(:));
orthoviewRGB(resimg);

figure(3);
result = (result>0).*result;
resimg(:,:,:,1) = result/max(result(:));
resimg(:,:,:,2) = labelimg>0.5; %/max(labelimg(:));
resimg(:,:,:,3) = 0*imgorig.*(imgorig>0) / max(imgorig(:));
orthoviewRGB(resimg);



model.L = L;
model.gwidin = gwidin;
model.gwidout = gwidout;
model.alpha = alpha;
model.ranktup = ranktup;
model.gamma = gamma;
model.inrank = inrank;


return


