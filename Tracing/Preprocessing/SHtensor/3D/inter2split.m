function split = inter2split(inter)
    sz = size(inter);
    sz(1) = sz(1)/2;
    N = length(inter(:));  
    split = inter(1:2:N) + 1i* inter(2:2:N);
    split = reshape(split,sz);
end