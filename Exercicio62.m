clear all; 
a=[-18.5:0.1:-17.5];
b=[-16.6:0.1:-16];
%b=[0.9:0.01:1.1];
 for x=1:length(b);
     for y=1:length(a);
         z(y,x)=fxy2([b(x) a(y)]);
     end;
 end;
 figure;
 mesh(z);
 figure;
 contour(z,50);