%{
***************************************************************************
* Otimização Convexa com 
* Steepest descent without line search
* conforme Antoniou Cap 5 - algoritmo 5.1 + Eq 5.17.
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
epsi=1e-6;
x0=[4 4]';
k = 1; xk = x0;
fk = fxy2(xk);
gk = dfxy2(xk);
ah = 1;
fh = fxy2(xk-ah*gk);
gk2 = gk'*gk;
ak = (gk2*ah^2)/(2*(fh-fk+ah*gk2));
adk = -ak*gk;
er = norm(adk);
while er >= epsi,
   xk = xk + adk;
   fk = fxy2(xk);
   gk = dfxy2(xk);
   ah = ak;
   fh = fxy2(xk-ah*gk);
   gk2 = gk'*gk;
   ak = (gk2*ah^2)/(2*(fh-fk+ah*gk2));
   adk = -ak*gk;
   er = norm(adk);
   k = k + 1;
end
xs = xk + adk; fs = fxy2(xs);
disp(sprintf('minimo=%2.8f em %d passos com precisão de %1.1e',fs,k,epsi));
disp(sprintf('x_otimo='));
xs