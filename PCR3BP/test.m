X = [0.8170501798 0.25 -0.38 -0.1409438404],

x = X(1); y = X(2); vx = X(3); vy = X(4); 

r1 = ((x+const.mu)^2+y^2)^(3/2);
r2 = ((x-1+const.mu)^2+y^2)^(3/2);

ax = 2*vy+x-((1-const.mu)*(x+const.mu)/r1)-((x-1+const.mu)*(const.mu)/r2);
ay = -2*vx+y-((1-const.mu)*y/r1)-((const.mu)*y/r2);

f = [vx vy ax ay],
