%{
***************************************************************************
* Otimização Convexa com 
* CVX - exercício 1 lista 4 de Otimização Convexa
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
A1=-1/8*[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; B1=[1; 2; 1; 2];
A2=-[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; B2=[1; 1; 2; 4];
warning off;
cvx_begin
    cvx_precision high
    variable x(4) nonnegative
    minimize (x(1)+1.5*x(2)+x(3)+x(4));
    subject to
        %-1/4*x(1)^2-x(2)^2+1/2*x(1)+3/4 >= 0;
        %(-5/8*x(3)^2)-(5/8*x(4)^2)-(6/8*x(3)*x(4))+(11/2*x(3))+(13/2*x(4))-(35/2) >= 0;
        x'*A1*x + x'*B1 - 3 == 0;
        x'*A2*x + x'*B2 - 5 == 0;
cvx_end
%minimize (x(1)^2 - 2*x(3)*x(1) + x(3)^2 + x(2)^2 -2*x(2)*x(4) + x(4)^2);
warning on;
format long
x