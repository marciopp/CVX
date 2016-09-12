%{
***************************************************************************
* Otimização Convexa
* Jacobiano
* conforme Antoniou Cap 5
* Marcio Pinto Pereira - julho de 2016
* Programado em Matlab R2016a 
* Licenciado sob CC-BY-SA
***************************************************************************
%}
function[y] = jdfxy2(x)
%y= -25*(x^4) +16*(x^3) -36*(x^2) +22*(x) -2;
%y= (-0.2/(x^0.8))-2*(log(10-x)/(10-x))+2*(log(x-2)/(x-2));
%y= [ 2.8*x(1)^3-x(2)*sin(x(1)*x(2))-16*x(1)-8; 12*x(2)-x(1)*sin(x(1)*x(2))];

%dfx=2*(2*x(1)*(x(1)^2+x(2)^2-1)+x(1)+x(2)-1);
%dfy=2*(2*x(2)*(x(1)^2+x(2)^2-1)+x(1)+x(2)-1);
%y=[dfx; dfy];
%{
dfx=2*(200*(x(1)-x(4))^3+x(1)+10*x(2));
dfy=4*(5*(x(1)+10*x(2))+(x(2)-2*x(3))^3);
dfz=10*(x(3)-x(4))-8*(x(2)-2*x(3))^3;
dfw=10*(40*(x(4)-x(1))^3+x(4)-x(3));

y=[dfx dfy dfz dfw];
%}

df1x=2*x(1)+20*x(2);
df2x=0;
df3x=0;
df4x=400*(x(1)-x(4))^3;
df1y=20*(10*x(2)+x(1));
df2y=0;
df3y=4*(x(2)-2*x(3))^3;
df4y=0;
df1z=0;
df2z=10*(x(3)-x(4));
df3z=8*(2*x(3)-x(2))^3;
df4z=0;
df1w=0;
df2w=10*(x(4)-x(3));
df3w=0;
df4w=400*(x(4)-x(1))^3;

y=[ df1x df1y df1z df1w ; df2x df2y df2z df2w ; df3x df3y df3z df3w ; df4x df4y df4z df4w ];

end