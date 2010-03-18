function Y = split2inter(re,im)
N = length(re(:));  
odd = 1:2:(2*N-1);
Y(odd) = re(:);
Y(odd+1) = im(:);
sz = size(re);
sz(1) = sz(1)*2;
Y = reshape(Y,sz);