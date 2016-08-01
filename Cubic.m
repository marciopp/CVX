%{
***************************************************************************
* Otimização Convexa com 
* Cubic Search
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
%% Cubic Search
cont=0; x0=10e99; xbarra=0; x1=xl; x3=xu; x2=(xu+xl)/2;
fx1=fx(xl); fx2=fx(x2);fx3=fx(x3);
f1linha=dfx(x1);
while (1)
    cont=cont+1; disp(cont);
    beta = (fx(x2)-fx(x1)+f1linha*(x1-x2))/((x1-x2)^2);
    gamma= (fx(x3)-fx(x1)+f1linha*(x1-x3))/((x1-x3)^2);
    theta= (2*x1^2-x2*(x1+x2))/(x1-x2);
    phi  = (2*x1^2-x3*(x1+x3))/(x1-x3);
    a3   = (beta-gamma)/(theta-phi);
    a2   = beta-theta*a3;
    a1   = f1linha-2*a2*x1-3*a3*x1^2;
    aux  = sqrt(a2^2-3*a1*a3);
    extremo1 = (-a2+aux)/(3*a3); disp(extremo1);
    extremo2 = (-a2-aux)/(3*a3); disp(extremo2);
    if (extremo1 > (-a2/(3*a3)))
        xbarra = extremo1;
    else
        xbarra = extremo2;
    end;
    fbarra=fx(xbarra);
    if (abs(xbarra-x0) < tolerancia)
        disp(sprintf('x_otimo=%2.8f, minimo=%2.8f em %d passos com precisão de %1.1e',xbarra,fx(xbarra),cont,tolerancia));
        return;
    end;
    [fmax, mfmax]= max([fx(x1) fx(x2) fx(x3)]);
    x0=xbarra; xm=xbarra; fmax=fbarra;
    if (mfmax == 1) f1linha=dfx(xbarra); end;
    
end;
