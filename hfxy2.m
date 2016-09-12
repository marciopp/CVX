%{
***************************************************************************
* Otimização Convexa
* Hessianas
* conforme Antoniou Caps 4 e 5
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
function[y] = hfxy2(x)
%y= -25*(x^4) +16*(x^3) -36*(x^2) +22*(x) -2;
%y= (-0.2/(x^0.8))-2*(log(10-x)/(10-x))+2*(log(x-2)/(x-2));
%y= [ 2.8*x(1)^3-x(2)*sin(x(1)*x(2))-16*x(1)-8; 12*x(2)-x(1)*sin(x(1)*x(2))]

%{
dfxx=12*x(1)^2+4*x(2)^2-2;
dfxy=8*x(1)*x(2)+2;
dfyx=8*x(1)*x(2)+2;
dfyy=12*x(2)^2+4*x(1)^2-2;
y=[dfxx dfxy; dfyx dfyy];
%}

dfxx=2*(600*(x(1)-x(4))^2+1);
dfxy=20;
dfxz=0;
dfxw=-1200*(x(1)-x(4))^2;
dfyx=20;
dfyy=4*(3*(x(2)-2*x(3))^2+50);
dfyz=-24*(x(2)-2*x(3))^2;
dfyw=0;
dfzx=0;
dfzy=-24*(x(2)-2*x(3))^2;
dfzz=48*(x(2)-2*x(3))^2+10;
dfzw=-10;
dfwx=-1200*(x(1)-x(4))^2;
dfwy=0;
dfwz=-10;
dfww=10*(120*(x(4)-x(1))^2+1);
y=[dfxx dfxy dfxz dfxw; dfyx dfyy dfyz dfyw; dfzx dfzy dfzz dfzw; dfwx dfwy dfwz dfww];
%}
end