function z=f(x,y)
if x<1 & x>-1 &y<1 & y>-1 
    z =x*x+y*y;
elseif (x==1|x==-1) & y<1&y>-1
    z=-y*y;
elseif x<1&x>-1&(y==1|y==-1)
    z=-x^2;
else
    z=0;
end
