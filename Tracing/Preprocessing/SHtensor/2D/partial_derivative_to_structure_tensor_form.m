function Sdata = partial_derivative_to_structure_tensor_form(I)

% partial_derivative_to_structure_tensor_form - sets up structure tensor form from Ix, Iy, ...%%%%
% 
% Example:
%   pDerData = [Ix, Iy, It];
%   S = partial_derivative_to_structure_tensor_form(pDerData);   %%% sets up the 3D structure tensor
% 
% Author: Shawn Arseneau
% Created: September 20, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = numel(I);
    Sdata = zeros(N, N);

    for e=1:N
        Sdata(e,e) = I(e)*I(e);            
    end

