 a=[-pi:0.1:pi];
 for x=1:length(a);
     for y=1:length(a);
         z(x,y)=0.7*(a(x)^4) -8*(a(x)^2) +6*(a(y)^2) -8*a(x) +cos(a(x)*a(y));
     end;
 end;
 mesh(z);