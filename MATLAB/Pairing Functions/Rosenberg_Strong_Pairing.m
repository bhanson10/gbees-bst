state = [100 200],

x = state(1); y = state(2);

key = (max(x,y))^2+max(x,y)+x-y;
key,

z = key; m = floor(key^(1/2));
if(z-m^2<m)
    x=z-m^2; y=m;
else
    x=m; y=m^2+2*m-z;
end

state = [x y],