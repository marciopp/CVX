%{
***************************************************************************
* Otimização Convexa com 
* Fletcher's Inexact Line Search
* conforme Antoniou Cap 4
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
clear all;
close all;
x0=[-pi pi]'; d0 =[1.0 -1.3]';
tau = 0.1; chi = 0.75; rho = 0.1; sigma = 0.1; epsilon = 1e-10;
f0=fxy2(x0); gk=dfxy2(x0); 
aL = 0; alphau = 1e99; deltaf0 = f0;
dk = d0; fL = f0; flinhaL = gk'*dk;
alpha0 = 1;
while 1,
    deltak = alpha0*dk;
    f0=fxy2(x0+deltak);
    if ((f0 > (fL + rho*(alpha0 - aL)*flinhaL)) & (abs(fL - f0) > epsilon))
      if (alpha0 < alphau)
        alphau = alpha0;
      end     
      alpha0til = aL + ((alpha0 - aL)^2*flinhaL)/(2*(fL - f0 + (alpha0 - aL)*flinhaL)); % Eq 4.57
      alpha0Ltil = aL + tau*(alphau - aL);
      if (alpha0til < alpha0Ltil)
        alpha0til = alpha0Ltil;
      end
      alpha0Util = alphau - tau*(alphau - aL);
      if (alpha0til > alpha0Util)
        alpha0til = alpha0Util;
      end
      alpha0 = alpha0til;
    else
      gtemp=dfxy2(x0+alpha0*dk);
      flinha0 = gtemp'*dk;
      if (((flinha0 < sigma*flinhaL) & (abs(fL - f0) > epsilon) & (flinhaL ~= flinha0) ))
        deltaalpha0 = (alpha0 - aL)*flinha0/(flinhaL - flinha0); % Eq 4.58
        if (deltaalpha0 <= 0)
          alpha0til = 2*alpha0;
        else
          alpha0til = alpha0 + deltaalpha0;
        end
        alpha0Util = alpha0 + chi*(alphau - alpha0);
        if (alpha0til > alpha0Util)
          alpha0til = alpha0Util;
        end
        aL = alpha0;
        alpha0 = alpha0til;
        fL = f0;
        flinhaL = flinha0;
      else
        break;
      end
    end
end
alpha0
 

 