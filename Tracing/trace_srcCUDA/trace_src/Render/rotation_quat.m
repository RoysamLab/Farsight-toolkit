function R = rotation_quat( Q )
% Gives the orthogonal matrix R corresponding to a rotation by a unit
% quaternion Q
w = Q(1);
x= Q(2);
y= Q(3);
z= Q(4);

%w = w/N;
%x = x/N;
%y = y/N;
%z = z/N;

R = [[1-2*(y^2+z^2) 2*(x*y-w*z) 2*(x*z+w*y)];
     [2*(x*y+w*z) 1-2*(x^2+z^2) 2*(y*z-w*x)];
     [2*(x*z-w*y) 2*(y*z+w*x) 1-2*(x^2+y^2)]];
