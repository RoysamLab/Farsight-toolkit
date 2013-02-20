function renderSElocal(s,N)
if nargin == 1
    N = 39;
end
t = SphericalTesselation(N);
t = GetSurface(t,s);
f = [t.n1' t.n2' t.n3'] ;
v = t.eb';
patch('Faces',f,'Vertices',v, 'facecolor','w');