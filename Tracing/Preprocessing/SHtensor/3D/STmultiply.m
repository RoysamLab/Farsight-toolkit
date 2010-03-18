
function STres = STmult(ST1,ST2,J,factor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  result = STmult(ST1,ST2,J,factor)
%
%  ST1,ST2 - Spherical tensor images
%  J       - Spherical rank of output (obeying trangle inequality)
%  factor  - multiply result with factor (optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4,
    factor = 1;
end;

L1 = (size(ST1,1)/2-1);
L2 = (size(ST2,1)/2-1);

if J > L1+L2 | J < abs(L1-L2),
    disp('triangle condition violated!');
    return;
end;


STres = STmult(ST1,ST2,J,myREAL(factor));
