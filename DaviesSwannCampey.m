%{
***************************************************************************
* Otimização Convexa com 
* Davies, Swann and Campey Search
* conforme Antoniou Cap 4
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
clear all;
close all;
% f(x)= -5x^5+4x^4-12x^3+11x^2-2x+1, entre -0.5 e 0.5
xl=-0.5; xu=0.5; % intervalo inicial 
delta1=0.1*xl; % 0.1 do ponto inicial
K=0.1; % K entre 0 e 1
p=1; % 1 ou -1
tolerancia=1e-6; % critério de parada
% f(x)= (ln(x-2))^2+(ln(10-x))^2-x^0.2, entre 6 e 9.9
% f(x)= -3x*sin(0.75x)+e^-2x, entre 0 e 2pi
% range of uncertainty = xl, xu
% epsilon: incerteza, x1 no meio entre xa e xb
% inicio em xa= x1 + epsilon/2 : xl < xa < xu
% e         xb= x1 - epsilon/2 : xl < xb < xu
% caso a: f(xa) < f(xb) : xl = x1 + epsilon/2
% caso b: f(xa) > f(xb) : xu = x1 - epsilon/2
% caso c: f(xa) = f(xb) : qq acima
% 
%% Davies, Swann and Campey Search
%fxa=fx(x1-dk); fxb=(x1+dk);
cont=0; x0=10e99;xbarra=0; x1=(xu+xl)/2; 
while (abs(xbarra-x0) > tolerancia)
    x0=xbarra;
    fxa=fx(xl); fxb=fx(x1);fxc=fx(xu);
    xaux1=(x1^2-xu^2)*fxa + (xu^2-xl^2)*fxb + (xl^2-x1^2)*fxc;
    xaux2=2*((x1-xu)*fxa + (xu-xl)*fxb + (xl-x1)*fxc);
    xbarra=xaux1/xaux2; fbarra=fx(xbarra);
    if ( xbarra < x1)
        if (fbarra <= fxb)
            xu = x1; fxc = fxa; x1 = xbarra; fxa = fbarra;
        else
            xl = xbarra; fxa = fbarra;
        end;
    else
        if (fbarra <= fxb)
            xl = x1; fxa = fxb; x1 = xbarra; fxa = fbarra;
        else
            xu = xbarra; fxc = fbarra;
        end;
    end;
    cont=cont+1; %disp(cont);
end;
disp(sprintf('x_otimo=%2.8f, minimo=%2.8f em %d passos com precisão de %1.1e',x1,fx(x1),cont,tolerancia));