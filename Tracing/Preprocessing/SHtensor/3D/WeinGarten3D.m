function [Tx, Ty,Tz, mag] = WeinGarten3D(vol)
%make sure that vol is smoothed before
W = smooth3(vol,'gaussian',5,2);

[d1, d2,d3] = size(vol);
[Gx,Gy,Gz] = gradient(vol);

Tx = zeros(d1,d2,d3);
Ty = zeros(d1,d2,d3);
Tz = zeros(d1,d2,d3);
mag = zeros(d1,d2,d3);
[Gxx Gyx Gzx] = gradient(Gx);
[Gxy Gyy Gzy] = gradient(Gy);
[Gxz Gyz Gzz] = gradient(Gz);
for i=2:d1-1
    i
    for j=2:d2-1
        for k=2:d3-1
            gx = Gx(i,j,k);
            gy = Gy(i,j,k);
            gz = Gz(i,j,k);

            gxx = Gxx(i,j,k);
            gxy = Gxy(i,j,k);
            gxz = Gxz(i,j,k);
            
            gyx = Gyx(i,j,k);
            gyy = Gyy(i,j,k);
            gyz = Gyz(i,j,k);

            gzx = Gzx(i,j,k);
            gzy = Gzy(i,j,k);
            gzz = Gzz(i,j,k);
            
            F2 = [gxx gxy gxz; 
                gyx gyy gyz; 
                gzx gzy gzz];
            F1 = [1+gx*gx , gx*gy , gx*gz; 
                gy*gx , 1+gy*gy , gy*gz; 
                gz*gx , gz*gy , 1+gz*gz];
            
            d = 1 + gx*gx + gy*gy + gz*gz;
            
            W = F2 * inv(F1) / sqrt(d);
            
            [u,v] = eig(W); 
            %assuming that you are taking minimum absolute ev
            [m,n] = min(abs(diag(v)));
           
            Tx(i,j,k) = u(1,n);
            Ty(i,j,k) = u(2,n);
            Tz(i,j,k) = u(3,n);
            mag(i,j,k) = -sum(v(:));
        end
    end
end
mag(mag < 0) = 0;

            