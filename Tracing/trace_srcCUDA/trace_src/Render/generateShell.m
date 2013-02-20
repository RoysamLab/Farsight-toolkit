function W = generateShell(Segment, NN)
if nargin == 1
    NN = 17;
end
m = 2*pi*(0:NN)./(NN);
%x-y ellipse
x1 = sin(pi/2)*cos(m);
y1 = sin(pi/2)*sin(m);
z1 = cos(pi/2)*ones(size(y1));

%x-z ellipse
x2 = sin(m)*cos(0);
y2 = sin(m)*sin(0);
z2 = cos(m);

%y-z ellipse
x3 = sin(m)*cos(pi/2);
y3 = sin(m)*sin(pi/2);
z3 = cos(m);

W = [x1 x2 x3; y1 y2 y3; z1 z2 z3];
% W = [x3; y3; z3];

rho = calc_rho(W(1,:), W(2,:), W(3,:), Segment.e1, Segment.e2);

W(1,:) = W(1,:)*Segment.a1;
W(2,:) = W(2,:)*Segment.a2;
W(3,:) = W(3,:)*Segment.a3;

W = ([1 1 1]'*rho).*W;

R = rotation_quat(Segment.q);
W = R*W + repmat(Segment.mu,1,size(W,2));

