function g = GaussCache(sz,wid,fftwplanner)

global myGaussCache

for k = 1:length(myGaussCache),
    if (min(myGaussCache(k).sz == sz) == 1) && (myGaussCache(k).wid == wid),
        g = myGaussCache(k).gauss;
        return;
    end;
end

display('gauss cache: building gauss')

[X Y Z] = ndgrid(0:(sz(1)-1),0:(sz(2)-1),0:(sz(3)-1));
X = X - ceil(sz(1)/2);
Y = Y - ceil(sz(2)/2);
Z = Z - ceil(sz(3)/2);
R2 = myREAL(X.^2 + Y.^2 + Z.^2);    
g = myREAL(fftshift(exp(-R2/wid^2)));

display('gauss cache: ffting gauss')

g = FFTfilterkernel(g,fftwplanner);

    
myGaussCache(end+1).sz = sz;
myGaussCache(end).wid = wid;
myGaussCache(end).gauss = g;



function gauss = FFTfilterkernel(gauss,fftwplanner)

    gauss = reshape(gauss,[1 size(gauss)]);
    gauss(2,:,:,:) = 0;
    gauss = myfftw(gauss,[2 3 4],fftwplanner);
    gauss(2,:,:,:) = gauss(1,:,:,:);
    gauss = gauss / gauss(1);


