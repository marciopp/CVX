clear all; close all; 
a=[-.5:0.01:1.2];
 for x=1:length(a);
     for y=1:length(a);
         z(x,y)=fxy2([a(y); a(x)]);
     end;
 end;
 figure;
 meshc(z);
 figure;
 contour(z,100);