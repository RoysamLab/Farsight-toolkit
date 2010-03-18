function out = mySTmultRealAsym(A,B,l)
    A = split2inter(real(A),imag(A));
    B = split2inter(real(B),imag(B));
    j1 = size(A,1)/2-1;
    j2 = size(B,1)/2-1;
    out = STmultiply(A,B,l);
    out = inter2split(out);
%       for co = 1:numel(B)
%           out(1,:,:,:) = B(co).*A;
%       end
end