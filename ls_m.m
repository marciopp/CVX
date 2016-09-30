function a = ls_m(fcname,xk,yk,dx,dy,ak,tau,bt)
aa = (0:0.005:1)*ak;
ps = zeros(101,1);
for i = 1:101,
    ai = aa(i);
    xki = xk + ai*dx;
    yki = yk + ai*dy;
    fci = feval(fcname,xki);
    fi = fci(1);
    cki = fci(2:end);
    yci = yki - cki;
    ps(i) = fi - tau*sum(log(yki)) + 0.5*bt*(yci'*yci);
end
[psm,ind] = min(ps);
a = aa(ind);
a = min(1,a);