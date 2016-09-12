%{
***************************************************************************
* Otimização Convexa com 
* Steepest descent
* conforme Antoniou Cap 5 - algoritmo 5.1 + Eq 5.17.
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
epsi=1e-6;
k = 1;
xk = [4 4]';
gk = dfxy2(xk);
dk = -gk;
ak = FILS(xk,dk);
adk = ak*dk;
er = norm(adk);
while er >= epsi,
   xk = xk + adk;
   gk = dfxy2(xk);
   dk = -gk;
   ak = FILS(xk,dk);
   adk = ak*dk;
   er = norm(adk);
   k = k + 1;
end
xs = xk + adk; fs = fxy2(xs);
disp(sprintf('minimo=%2.8f em %d passos com precisão de %1.1e',fs,k,epsi));
disp(sprintf('x_otimo='));
xs
