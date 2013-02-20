function p = calc_rho(x,y,z,e1,e2)

%p=((x.^2+y.^2).^(1/e1)+abs(z).^(2/e1)).^(-(e1/2));
p=exp((-(e1/2))*log(exp(1/e1*log(x.*x+y.*y))+exp((1/e1)*log(z.*z))));
%p=((abs(z).^(2/e2)+abs(y).^(2/e2)).^(e2/e1)+abs(x).^(2/e1)).^(-(e1/2));
