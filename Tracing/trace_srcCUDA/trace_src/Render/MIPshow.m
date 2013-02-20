function h = MIPshow(Vol, F)
if nargin == 1
    F = 'min';
end
[d1,d2,d3] = size(Vol);
Myx = feval(F,Vol,[],3);
Mxz = feval(F,Vol,[],1);
Mxz = reshape(Mxz,d2,d3);
Myz = feval(F,Vol,[],2);
Myz = reshape(Myz,d1,d3);
% M = [imrotate(mat2gray(Myz),90) zeros(d3,d3); mat2gray(Mxy) (mat2gray(Mxz))];
M = [Myx, Myz ; Mxz', zeros(d3,d3) ];
% Vol = Vol(Vol>0);
%Vol = Vol(Vol<255);
maxVol = max(Vol(:));
minVol = min(Vol(:));
figure(1);
imagesc(M, [minVol, maxVol]);
% set(gca,'Ydir','normal');
colormap(gray);
axis image;
colorbar;
title(['Max value: ', num2str(maxVol), 'Min value: ', num2str(minVol)]);
ylabel('Y (Dim 1)                      Z (Dim 3)');
xlabel('X (Dim 2)                      Z (Dim 3)');


