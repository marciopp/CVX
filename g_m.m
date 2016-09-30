function z = g_m(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
Ak = [eye(4); -eye(4)];
dx1=400*x1^3-200*x1*x2+2*x1-2;
dx2=-100*x1^2+200*x2+20.2*x2-20.2+19.8*x4-19.8;
dx3=360*x3^3-180*x3*x4+2*x3-2;
dx4=-90*x3^2+180*x4+20.2*x4-20.2+19.8*x2-19.8;
gk = [dx1 ; dx2 ; dx3 ; dx4];
z = [gk Ak'];