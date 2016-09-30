%{
***************************************************************************
* Otimização Convexa com 
* conforme Antoniou - cap 12 - alg 12.6
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* baseado no livro-texto e códigos 
* de Practical Optimization - Antoniou e Lu
* Licenciado sob CC-BY-SA
***************************************************************************
%}
c=[1 1.5 1 1]';
A=[1 2 1 2; 1 1 2 4];
b=[3; 5];
x0=[1 1 1 1]';
lmd0=[1 1]';
mu0=[1 1 1 1]';
epsi=1e-6;
n = length(x0);
x = x0(:);
lam = lmd0(:);
mu = mu0(:);
gap = mu'*x;
k = 0;
while gap > epsi,
  tau_h = gap/n;
  mui = 1./mu;
  D = diag(x./mu);
  Y = inv(A*D*A');
  rd = c - A'*lam - mu;
  d_lam_a = Y*(b+A*D*rd);
  d_mu_a = rd - A'*d_lam_a;
  d_x_a = -x - D*d_mu_a;
  a_p = [];
  a_d = [];
  for i = 1:n,
    if d_x_a(i) < 0,
     a_p = [a_p x(i)/(-d_x_a(i))];
    end
    if d_mu_a(i) < 0,
     a_d = [a_d mu(i)/(-d_mu_a(i))];
    end
  end
  a_p_aff = min([1 min(a_p)]);
  a_d_aff = min([1 min(a_d)]);
  tau_aff = ((mu + a_d_aff*d_mu_a)'*(x + a_p_aff*d_x_a))/n;
  sigma_k = (tau_aff/tau_h)^3;
  tau = sigma_k*tau_h;
  y = mui.*d_x_a.*d_mu_a - tau*mui;
  d_lam_c = Y*A*y;
  d_mu_c = -A'*d_lam_c;
  d_x_c = -y - D*d_mu_c;
  d_x = d_x_a + d_x_c;
  d_mu = d_mu_a + d_mu_c;
  d_lamda = d_lam_a + d_lam_c;
  a_p = [];
  a_d = [];
  for i = 1:n,
    if d_x(i) < 0,
     a_p = [a_p x(i)/(-d_x(i))];
    end
    if d_mu(i) < 0,
     a_d = [a_d mu(i)/(-d_mu(i))];
    end
  end
  a_k_p = min([1 0.99*min(a_p)]);
  a_k_d = min([1 0.99*min(a_d)]);
  x = x + a_k_p*d_x;
  lam = lam + a_k_d*d_lamda;
  mu = mu + a_k_d*d_mu;
  gap = mu'*x;
  k = k + 1;
end
xs = x
fs = c'*x