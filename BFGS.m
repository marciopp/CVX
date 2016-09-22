%{
***************************************************************************
* Otimização Convexa com 
* BFGS
* conforme Antoniou - cap 6
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
epsi = 1e-6;
x0 = [0 0]';
n = length(x0);
I = eye(n);
k = 1;
xk = x0;
Sk = I;
fk = fxy2(xk);
gk = dfxy2(xk);
dk = -Sk*gk;
ak = FILS(xk,dk);
dtk = ak*dk;
xk_new = xk + dtk;
fk_new = fxy2(xk_new);
dfk = abs(fk - fk_new);
err = max(dfk,norm(dtk));
while err >= epsi,
      gk_new = dfxy2(xk_new);
      gmk = gk_new - gk;
      D = dtk'*gmk;
      if D <= 0,
         Sk = I;
      else
         sg = Sk*gmk;
         sw0 = (1+(gmk'*sg)/D)/D;
         sw1 = dtk*dtk';
         sw2 = sg*dtk';
         Sk = Sk + sw0*sw1 - (sw2'+sw2)/D;
      end
      fk = fk_new;
      gk = gk_new;
      xk = xk_new;
      dk = -Sk*gk;
      ak = FILS(xk,dk);
      dtk = ak*dk;
      xk_new = xk + dtk;
      fk_new = fxy2(xk_new);
      dfk = abs(fk - fk_new);
      err = max(dfk,norm(dtk));
      k = k + 1;
end
xs = xk_new
