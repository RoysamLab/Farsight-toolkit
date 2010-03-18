function [dir,J] = get_weingarten_in_spherical(im)

epsi = 0.0000001;
[Vx,Vy,Vz,J] = WeinGarten3D(single(im));
maxJ = max(J(:));
J = J./maxJ;
h = vol3d('cdata',J,'texture','3D');
vol3d(h);
dir = zeros([2 size(im)]);
sign1 = sign(rand(size(Vx))-0.5);
sign1 = sign(Vx);
dir(1,:,:,:) = - 10*Vx.*sign1.*J/sqrt(2)+epsi - 1i*sign1.*(10*Vy.*J/sqrt(2)+epsi) ;
dir(2,:,:,:) = 10*sign1.*Vz.*J+epsi;
dir = split2inter(real(dir),imag(dir));

end