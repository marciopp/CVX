%{
***************************************************************************
* Otimização Convexa
* Funções estudadas
* conforme Antoniou Cap 4 e 5
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
function[y] = fxy2(x)
%y= -5*(x^5) +4*(x^4) -12*(x^3) +11*(x^2) -2*x +1;
%y= (log(x-2))^2+(log(10-x))^2-x^0.2;
%y= -3*x*sin(0.75*x)+exp(-2*x);
%y=0.7*x(1)^4-8*x(1)^2+6*x(2)^2+cos(x(1)*x(2))-8*x(1);
%y=(x(1)^2+x(2)^2-1)^2+(x(1)+x(2)-1)^2;
y=(x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+100*(x(1)-x(4))^4;
end