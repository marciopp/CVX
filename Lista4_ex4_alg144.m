%{
***************************************************************************
* Otimização Convexa com 
* conforme Antoniou - cap 14 - alg 14.4
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* baseado no livro-texto e códigos 
* de Practical Optimization - Antoniou e Lu
* Licenciado sob CC-BY-SA
***************************************************************************
%}
c = [1 0 2 -1]';
[F0,FF] = lista4_4_FF;
epsi = 1e-6;
[m,nn] = size(FF);
p = round(nn/m);
pz = p + 1;
mz = m + 1;
ct = [c(:); 0];
dt = [zeros(p,1); 1];
FFe = [FF F0];
FFt = zeros(mz,pz*mz);
FFt(mz,end) = 1;
for i = 1:pz,
    FFt(1:m,(i-1)*mz+1:i*mz-1) = FFe(:,(i-1)*m+1:i*m);
end
I = eye(mz);
Xki = I;
[x,Xkp] = proj(Xki,FFt,mz,pz);
vi = min(eig(Xkp));
while vi <= 0,
   Xw = Xki*Xkp - I;
   rou = max(abs(eig(Xw)));
   gk = 1/(1+rou);
   Xki = Xki - gk*Xw*Xki;
   [x,Xkp] = proj(Xki,FFt,mz,pz);
   vi = min(eig(Xkp));
end
k = 0;
Xk = Xkp;
xk = x;
fk = (ct'*xk)/(dt'*xk);
Xki = inv(Xk);
Xkp = Xk;
ve = min(eig(Xkp));
er = 1;
while er >= epsi,
   Wk = Xki*Xkp - I;
   if ve <= 0,
      Yk = Xkp - Xk;
   else
      rkd = trace(Wk*Wk) - 0.99^2;
      if rkd >= 0,
         Yk = Xkp - Xk;
      else
         [Fx,Fwi,xc,xd] = find_cd(ct,Xki,FFt,mz,pz);
         c2 = ct'*xc;
         d2 = dt'*xd;
         e2 = ct'*xd;
         ck = ct'*xk;
         dk = dt'*xk;
         wa = dk^2 + rkd*d2;
         wb = -(ck*dk+e2*rkd);
         wc = ck^2 + rkd*c2;
         fw = sqrt(wb^2 - wa*wc);
         fk = (-wb - fw)/wa;
         qw = zeros(pz,1);
         for i = 1:pz,
             qw(i) = trace(Fx(:,(i-1)*mz+1:i*mz));
         end
         pw = ct - fk*dt;
         pf = pw'*Fwi;
         qf = Fwi*qw;
         lams = -(pf*qw)/(pf*pw);
         xs = lams*pf' + qf;
         Xkpf = zeros(mz);
         for i = 1:pz,
             Xkpf = Xkpf + xs(i)*FFt(:,(i-1)*mz+1:i*mz);
         end
         Yk = Xkpf - Xk;
      end
   end
   rou = max(abs(eig(Xki*Yk)));
   gk = 1/(1+rou);
   Xki = Xki - gk*Xki*Yk*Xki;
   Xk = inv(Xki);
   [xk,Xkp] = proj(Xki,FFt,mz,pz);
   ve = min(eig(Xkp));
   vf = (ct'*xk)/(dt'*xk);
   k = k + 1;
   er = vf - fk;
end
xs = xk(1:p)/xk(p+1)
fs = c'*xs
Fs = Xkp(1:m,1:m)