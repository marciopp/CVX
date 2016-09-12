%{
***************************************************************************
* Otimização Convexa com 
* Gauss-Newton e Fletcher's Inexact Line Search
* conforme Antoniou Cap 5 - algoritmo 5.5
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
epsi=1e-6; x0=[200 -200 100 -100]';
%x0=[-2 -1 1 2]';
x = x0(:);
n = length(x);
In = eye(n);
k = 1;
xk = x0;
F_k = fxy2(xk);
gk = dfxy2(xk);
Jk = jdfxy2(xk);
Hk = 2*Jk'*Jk + 1e-12*In;

[L,D]=jfd4(Hk); % decomposição LDLt
yk=-L*gk;
dk = L'*inv(D)*yk;
%}
%dk = -inv(Hk)*gk;
ak = FILS(xk,dk);

adk = ak*dk;
er = norm(adk);
while er > epsi,
  xk = xk + adk;
  F_k1 = fxy2(xk);
  gk = dfxy2(xk);
  Jk = jdfxy2(xk);
  Hk = 2*Jk'*Jk + 1e-12*In;
    
  [L,D]=jfd4(Hk);
  yk=-L*gk;
  dk = L'*inv(D)*yk;
  %}
  %dk = -inv(Hk)*gk;
  ak = FILS(xk,dk);
  adk = ak*dk;
  er = abs(F_k1 - F_k);
  k = k + 1;
  F_k = F_k1;
end
xs = xk + adk; fs = fxy2(xs);
disp(sprintf('minimo=%2.8f em %d passos com precisão de %1.1e',fs,k,epsi));
disp(sprintf('x_otimo='));
xs