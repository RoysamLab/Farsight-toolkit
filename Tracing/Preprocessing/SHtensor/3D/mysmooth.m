

function out = mysmooth(inputimg,gauss,fftwplanner)
    
    imgFT = myfftw(inputimg,[2 3 4],fftwplanner);
    imgFT = imgFT .* gauss;  
    out = myifftw(imgFT,[2 3 4],fftwplanner);

