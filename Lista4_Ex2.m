%{
***************************************************************************
* Otimização Convexa com 
* CVX - exercício 2 lista 4 de Otimização Convexa
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
A1=-1/8*[0 0 0 0; 0 0 0 0; 0 0 5 3; 0 0 3 5]; B1=[0; 0; 11/2; 13/2];
A2=-[1/4 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]; B2=[1/2; 0; 0; 0];
warning off;
cvx_begin
    cvx_precision high
    variable x(4)
    minimize ((x(1)-x(3))^2 + (x(2)-x(4))^2);
    subject to
        %-1/4*x(1)^2-x(2)^2+1/2*x(1)+3/4 >= 0;
        %(-5/8*x(3)^2)-(5/8*x(4)^2)-(6/8*x(3)*x(4))+(11/2*x(3))+(13/2*x(4))-(35/2) >= 0;
        x'*A1*x + x'*B1 - 35/2 >= 0;
        x'*A2*x + x'*B2 + 3/4 >= 0;
cvx_end
%minimize (x(1)^2 - 2*x(3)*x(1) + x(3)^2 + x(2)^2 -2*x(2)*x(4) + x(4)^2);
warning on;
format long
x