%{
***************************************************************************
* Otimiza��o Convexa com 
* conforme Antoniou - cap 15 - alg 15.6
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* baseado no livro-texto e c�digos 
* de Practical Optimization - Antoniou e Lu
* Licenciado sob CC-BY-SA
***************************************************************************
%}
x0 = [1 1 1 1]';
%x0 = [0 0 0 0]';
%x0 = [2 2 2 2]';
%x0 = [10 10 10 10]';
y0 = [10 10 10 10 10 10 10 10]';
lmd0 = [1 1 1 1 1 1 1 1]';
tau0 = 0.001;
ei = 5e-4;
eo = 1e-5;
xk = x0(:);
yk = y0(:);
lmdk = lmd0(:);
tau = tau0;
n = length(x0);
q = length(y0);
q1 = q + 1;
L = 0;
K = 0;
do = 1;
E2 = eye(n);
while do >= eo,
   xk0 = xk;
   yk0 = yk;
   lmdk0 = lmdk;
   k = 0;
   di = 1;
     while di >= ei,
     bt = 0;
     fck = f_m(xk);
     ck = fck(2:q1);
     Gk = g_m(xk);
     gk = Gk(:,1);
     Ak = Gk(:,2:q1)';
     yki = 1./yk;
     wwk = yki.*lmdk;
     Yki = diag(yki);
     Lk = diag(lmdk);
     xlk = [n;q;xk;lmdk];
     Hk = h_m(xlk);
     Nk = Hk + Ak'*(Yki*Lk)*Ak;
     ew = min(eig(Nk));
     if ew <= 0,
        eb = -1.2*ew + 1e-6;
        Nk = Nk + eb*E2;
     end
     Nki = inv(Nk);
     gak = tau*yki - lmdk;
     rk = yk - ck;
     dx = Nki*(-gk + tau*Ak'*yki + Ak'*(wwk.*rk));
     dL = wwk.*(rk - Ak*dx) + gak;
     dy = Ak*dx - rk;
     ksi = gk - tau*Ak'*yki;
     sk = -ksi'*Nki*ksi + tau*yki'*rk + ksi'*(Nki*Ak')*(wwk.*rk);
     if sk >= 0,
        bmin  = sk/(norm(rk)^2);
        bt = 10*bmin;
     end
     awk = max([-dy./yk; -dL./lmdk]);
     ak = 0.95/awk;
     aks = ls_m('f_m',xk,yk,dx,dy,ak,tau,bt);
     xk = xk + aks*dx;
     yk = yk + aks*dy;
     lmdk = lmdk + aks*dL;
     k = k + 1;
     di = aks*(norm(dx) + norm(dy) + norm(dL));
  end
  tsi = q*min(yk.*lmdk)/(yk'*lmdk);
  tau = 0.1*(min([0.05*(1-tsi)/tsi,2]))^3*(yk'*lmdk)/q;
  K = k + K;
  L = L + 1;
  do = norm(xk-xk0) + norm(yk-yk0) + norm(lmdk-lmdk0);
end
xs = xk
fck = f_m(xk);
fs = fck(1)
k = K