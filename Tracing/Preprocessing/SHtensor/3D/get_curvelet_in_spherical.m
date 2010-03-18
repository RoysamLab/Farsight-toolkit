function [curveness] = get_curvelet_in_spherical(im)

epsi = 0.0000001;
[out,cos1,sin1] = get_curvelet(im);
out = double(out)./double(max(out(:)));
cos1 = double(cos1);
sin1 = double(sin1);

% h = vol3d('cdata',J,'texture','3D');
% vol3d(h);
curveness = zeros([2 size(im)]);

curveness(1,:,:,:) = - sin1.*out/sqrt(2) + 1i*out.*(cos1/sqrt(2)) ; % works in spiral
% curveness(1,:,:,:) =  cos1.*out/sqrt(2) + 1i*out.*(sin1/sqrt(2)) ; % does not work
curveness(2,:,:,:) = 0;
curveness = split2inter(real(curveness),imag(curveness));
% 
% imshow(squeeze(out(:,:,25)));hold on;
% [X,Y] = meshgrid(1:size(out,2),1:size(out,1));
% quiver(X,Y,squeeze(cos1(:,:,25)),squeeze(sin1(:,:,25)));



end