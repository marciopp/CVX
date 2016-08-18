%{
***************************************************************************
* Otimização Convexa com 
* Backtracking Line Search
* conforme Antoniou Cap 4
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
clear all;
close all;
% f(x)= -5x^5+4x^4-12x^3+11x^2-2x+1, entre -0.5 e 0.5
xl=0; xu=2*pi; % intervalo inicial 
alfa=0.15; beta=0.5; % alfa entre 0.01 e 0.3, beta entre 0.1 e 0.8
t=1;
tolerancia=1e-5; % critério de parada
% f(x)= (ln(x-2))^2+(ln(10-x))^2-x^0.2, entre 6 e 9.9
% f(x)= -3x*sin(0.75x)+e^-2x, entre 0 e 2pi
%% Backtracking Line Search
cont=0;
x=xl; % ponto inicial
while( t > tolerancia )
    %Determine a direção descendente
    delta_x=-sign(dfx(x)); % estamos em R
    while ((fx(x+t*delta_x)) > (fx(x)+alfa*t*dfx(x)*delta_x))
        t=beta*t;
        cont=cont+1; %disp(cont);
    end;
    x=x+t*delta_x;
end;
%end;
disp(sprintf('x_otimo=%2.8f, minimo=%2.8f em %d passos com precisão de %1.1e',x,fx(x),cont,tolerancia));