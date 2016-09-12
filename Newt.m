%{
***************************************************************************
* Otimização Convexa com 
* Newton com Fletcher's Inexact Line Search
* conforme Antoniou Cap 5 - algoritmo 5.3
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
epsi=1e-6; dt=.1;k = 1;
xk = [200 -200 100 -100]';
n = length(xk); In = eye(n,n);
gk = dfxy2(xk); Hk = hfxy2(xk); [V,D] = eig(Hk);
di = diag(D); dmin = min(di);
if dmin > 0,
   Hki = V*diag(1./di)*V';
else
   bt = dt - dmin;
   Hki = V*diag((1+bt)./(di+bt))*V';
end
dk = -Hki*gk;
ak = FILS(xk,dk);
adk = ak*dk;
er = norm(adk);
while er >= epsi,
   xk = xk + adk;
   gk = dfxy2(xk);
   Hk = hfxy2(xk);
   [V,D] = eig(Hk);
   di = diag(D);
   dmin = min(di);
   if dmin > 0,
      Hki = V*diag(1./di)*V';
   else
      bt = dt - dmin;
      Hki = V*diag((1+bt)./(di+bt))*V';
   end
   dk = -Hki*gk;
   ak = FILS(xk,dk);
   adk = ak*dk;
   er = norm(adk);
   k = k + 1;
end
xs = xk + adk; fs = fxy2(xs);
disp(sprintf('minimo=%2.8f em %d passos com precisão de %1.1e',fs,k,epsi));
disp(sprintf('x_otimo='));
xs