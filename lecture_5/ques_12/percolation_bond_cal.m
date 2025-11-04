syms p
f = p.^6 + 6.*p.^5.*(1-p) + 13.*p.^4.*(1-p).^2 + 10.*p.^3.*(1-p).^3 + 2.*p.^2.*(1-p).^4;
simplify(f)
g = simplify(diff(f,p));
gfun = @(p) 10.*p.^4 - 20.*p.^3 + 6.*p.^2 + 4.*p;
nu = log(2)/log(gfun(0.5));
disp(['临界值数计算',num2str(nu)])