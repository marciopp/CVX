%{
***************************************************************************
* Otimização Convexa com 
* Fletcher-Reeves
* conforme Antoniou - cap 6
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
epsi = 3e-7;
k = 1;
xk = [1 1]';
gk = dfxy2(xk);
dk = -gk;
g2 = gk'*gk;
ak = FILS(xk,dk);
adk = ak*dk;
err = norm(adk);
while  err >= epsi,
    xk = xk + adk;
    gk = dfxy2(xk);
    g2_new = gk'*gk;
    bk = g2_new/g2;
    dk = -gk + bk*dk;
    g2 = g2_new;
    ak = FILS(xk,dk);
    adk = ak*dk;
    err = norm(adk);
    k = k + 1;
end
xs = xk + adk
