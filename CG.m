%{
***************************************************************************
* Otimização Convexa com 
* Conjugate-gradient
* conforme Antoniou - cap 6
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
epsi=1e-6; k = 1;
xk = zeros(16,1);
gk = Exercicio61_g(xk);
Hk = Exercicio61_Q(xk);
dk = -gk;
g2 = gk'*gk;
ak = g2/(dk'*Hk*dk);
adk = ak*dk;
err = norm(adk);
while  err >= epsi,
    xk = xk + adk;
    gk = Exercicio61_g(xk);
    Hk = Exercicio61_Q(xk);
    g2_new = gk'*gk;
    bk = g2_new/g2;
    dk = -gk + bk*dk;
    g2 = g2_new;
    ak = g2/(dk'*Hk*dk);
    adk = ak*dk;
    err = norm(adk);
    k = k + 1;
end
xs = xk + adk
