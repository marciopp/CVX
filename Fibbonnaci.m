%{
***************************************************************************
* Otimização Convexa com 
* Fibonnacci Search
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
epsilon=1e-10; % quanto menor melhor. Limita a tolerância
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
%% Fibonaci - 1,1,2,3,5,8,13
cont=0; %fibonacci=zeros(100);
fn=(xu-xl)/tolerancia;
fibonacci(1)=1; fibonacci(2)=1; n=2;
while (fibonacci(n) < fn)
    n=n+1;
    fibonacci(n)=fibonacci(n-2)+fibonacci(n-1);
end;
for k=n:-1:3
    f=fibonacci(k-1)/fibonacci(k);
    xu2=xl+((xu-xl)*f); xl2=xu-((xu-xl)*f);
    fxu2=fx(xu2); fxl2=fx(xl2); % 2 cálculos pedidos ao oráculo
    if (fxu2 > fxl2)
        xu = xu2;
    else
        xl = xl2;
    end;
    cont=cont+1; %disp(cont);
end;
disp(sprintf('x_otimo=%2.8f, minimo=%2.8f em %d passos com precisão de %1.1e',(xu+xl)/2,fx((xu+xl)/2),cont,tolerancia));