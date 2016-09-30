%{
***************************************************************************
* Otimização Convexa com 
* conforme Antoniou - cap 13 - alg 13.3
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* baseado no livro-texto e códigos 
* de Practical Optimization - Antoniou e Lu
* Licenciado sob CC-BY-SA
***************************************************************************
%}
H = [4 0 0; 0 1 -1; 0 -1 1];
p = [-8 -6 -6]';
A = [1 1 1];
b = 3;
x0 = [1 2 2]';
lmd0 = -1;
mu0 = [0.2 0.2 0.2]';
rho = 3;
epsi = 1e-6;
n = length(x0);
rho_n = n + rho;
a_max = 1-1e-6;
x = x0(:);
lmd = lmd0(:);
mu = mu0(:);
xm = x.*mu;
gap = sum(xm);
k = 0;
while gap > epsi,
  tau = gap/rho_n;
  r_p = b - A*x;
  r_d = H*x + p - A'*lmd-mu;
  M = diag(mu);
  X = diag(x);
  Gam = inv(M+X*H);
  GaX = Gam*X*A';
  Y0 = inv(A*GaX);
  yd = Gam*(x.*(mu+r_d) - tau);
  d_lmd = Y0*(A*yd + r_p);
  d_x = GaX*d_lmd - yd;
  d_mu = H*d_x - A'*d_lmd + r_d;
  ind = find(d_x < 0);
  a_p = min(x(ind)./(-d_x(ind)));
  ind = find(d_mu < 0);
  a_d = min(mu(ind)./(-d_mu(ind)));
  a_k = a_max*min([a_p a_d]);
  x = x + a_k*d_x;
  lmd = lmd + a_k*d_lmd;
  mu = mu + a_k*d_mu;
  xm = x.*mu;
  gap = sum(xm);
  k = k + 1;
end
xs = x
fs = 0.5*xs'*(H*xs + 2*p)