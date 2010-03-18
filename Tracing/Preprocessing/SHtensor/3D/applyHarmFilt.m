function result = applyHarmFilt(imgorig,model,padsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  result = applyHarmFilt(img,model,padsize);
%
%  img     - the filter is applied on img
%  model   - a trained model
%  padsize - (optional) to pad img to padsize
%  result  - filter response
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fftwplanner = 1; %% fftwplan flag

%%%%%%%%%% padding
if nargin == 3,
    imgsz = padsize;
    unpadsz = size(imgorig);
    img = zeros(imgsz);
    img(1:unpadsz(1),1:unpadsz(2),1:unpadsz(3)) = imgorig;
else,
    imgsz = size(imgorig);
    unpadsz = size(imgorig);
    img = imgorig;
end;


gaussin = GaussCache(imgsz,model.gwidin,fftwplanner);
gaussout = GaussCache(imgsz,model.gwidout,fftwplanner);

fprintf('blurring input');
tic
inputimg = myREAL(img);
inputimg = reshape(inputimg,[1 imgsz]);
inputimg(2,:,:,:) = 0; 
inputimg_sav = inputimg;
inputimg = mysmooth(inputimg,gaussin,fftwplanner)/(prod(imgsz));
fprintf(' in %f s\n',toc);


display('computing stddev');
tic
simg = mysmooth(inputimg_sav.^2,gaussin,fftwplanner)/prod(imgsz);
stddev = simg(1,:,:,:) - inputimg(1,:,:,:).^2;
stddev = real(sqrt(stddev.* (stddev>0)));
toc




fprintf('computing derivatives');
tic
L = model.L;
D{1} = inputimg;
inrank(1) = 0;
for k = 2:L,
    D{k} = STderivUp(D{k-1});
    if k == 2 && model.gamma >= 0,          
            normD2 = model.gamma + stddev;
            D{2}(1,:,:,:) =  D{2}(1,:,:,:) ./normD2;
            D{2}(2,:,:,:) =  D{2}(2,:,:,:) ./normD2;
            D{2}(3,:,:,:) =  D{2}(3,:,:,:) ./normD2;
    end;   
    inrank(k) = k-1;
end;
%D{L+1}(1,:,:,:) = ones(imgsz);
%D{L+1}(2,:,:,:) = 0;
%inrank(L+1) = 0;
%D{1}(1,:,:,:) = 1; 
fprintf(' in %f s\n',toc);


% sort the output order of the products in ascending order
ranktup = model.ranktup;
Lmax = max(ranktup(:,3));
for l = 0:Lmax,
    ranks{l+1} = find(ranktup(:,3)==l)';
end;


display('applying filter');
result = 0;
for l = Lmax:-1:0,
    tic
    for j = ranks{l+1},
        result = result + STmultiply(D{ranktup(j,1)},D{ranktup(j,2)},ranktup(j,3),model.alpha(j));
        fprintf('.');
    end;
    if l >= 1,
        result = STderivDown(result);
    end;
    fprintf(' %i/%i  %f \n',Lmax-l+1,Lmax+1,toc);
end;



fprintf('smoothing output');
tic
result = mysmooth(result,gaussout,fftwplanner) /(prod(imgsz));
fprintf(' in %f s\n',toc);





result = reshape(result(1,1:unpadsz(1),1:unpadsz(2),1:unpadsz(3)),unpadsz);
result = (result>0).*result;



img2show = reshape(inputimg(1,:,:,:),size(result));
figure(2);
result = (result>0).*result; %% cut negative values
result = (result<1).*result + (result>=1); %% cut values greater 1

resimg(:,:,:,1) = result;
resimg(:,:,:,2) = img2show.*(img2show>0) / (max(img2show(:))+0.000001);
resimg(:,:,:,3) = 0;
orthoviewRGB(resimg);


figure(3)
orthoview(result/(max(result(:))+0.000001));

return

