%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
T=1; G.thresh=0.00002; G.max=10000; G.start=[-15; -4.5; 42]; 
G.dt=.00005; dt=.005; G.dx=0.5; G.d=3; G.sigma=10; G.b=(8/3); G.r=28;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
figure(1); clf; y=G.start; ys=y; 
for timestep=1:10000
  k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
  ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
end
plot3(ys(1,:),ys(2,:),ys(3,:),'g-','linewidth',1); view(-109,14);  hold on;
lighting phong; light('Position',[-1 0 0]); drawnow;

%{
figure(3); clf; y=G.start; ys=y; 
for timestep=1:10000
  k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
  ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
end
plot3(ys(1,:),ys(2,:),ys(3,:),'g-','linewidth',1); view(-109,14);  hold on;
lighting phong; light('Position',[-1 0 0]); drawnow;
%}

figure(2); clf; y=G.start; ys=y; 
for timestep=1:T/G.dt
  k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
  ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew;
end
view(-109,14);  hold on;
plot3(ys(1,:),ys(2,:),ys(3,:),'k-','linewidth',2); view(-109,14);  hold on;
plot3(ys(1,1),ys(2,1),ys(3,1),'k*'),plot3(ys(1,end),ys(2,end),ys(3,end),'k*');
lighting phong; light('Position',[-1 0 0]); drawnow;

for P=1:100; y=G.start+0.5*randn(3,1); ys=y;
  for timestep=1:T/G.dt
    k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
    ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew;
  end
  plot3(ys(1,:),ys(2,:),ys(3,:),'c-.','linewidth',0.3); 
  for i=1:round(T/(5*G.dt)):round(T/G.dt)+1
    plot3(ys(1,i),ys(2,i),ys(3,i),'k+');
    drawnow;
  end
end

%{
figure(4); clf; y=G.start; ys=y; 
for timestep=1:T/G.dt
  k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
  ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew;
end
view(-109,14);  hold on;
plot3(ys(1,:),ys(2,:),ys(3,:),'k-','linewidth',2); view(-109,14);  hold on;
plot3(ys(1,1),ys(2,1),ys(3,1),'k*'),plot3(ys(1,end),ys(2,end),ys(3,end),'k*');
lighting phong; light('Position',[-1 0 0]); drawnow;

for P=1:100; y=G.start+0.5*randn(3,1); ys=y;
  for timestep=1:T/G.dt
    k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
    ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew;
  end
  plot3(ys(1,:),ys(2,:),ys(3,:),'c-.','linewidth',0.3); 
  plot3(ys(1,1),ys(2,1),ys(3,1),'k+'), plot3(ys(1,end),ys(2,end),ys(3,end),'k+'), 
  plot3(ys(1,41),ys(2,41),ys(3,41),'k+'); plot3(ys(1,81),ys(2,81),ys(3,81),'k+');
  plot3(ys(1,121),ys(2,121),ys(3,121),'k+'); plot3(ys(1,161),ys(2,161),ys(3,161),'k+');
  drawnow;
end
%}

for j=0:0
    for i=0:5
        %file = "cpp_sim_data_v7_" + string(i) + ".txt"; 
        file = "pdf_"+string(j)+ "-" + string(i*round(T/(5*G.dt))) + ".txt";
        blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [4 inf]);
        D.P = blob(1,:); D.j = blob(2:4,:)'; D.n = length(D.P);
        Rotate_Plot(D,G,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(y,G)                          
    %f=[G.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -G.b*y(3)+y(1)*y(2)-G.b*G.r];
    f=[G.sigma*(y(2)-y(1));  y(1)*(G.r-y(3))-y(2);  y(1)*y(2)-G.b*y(3)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rotate_Plot(D,G,n)        
max_val = zeros(1,G.d); min_val = zeros(1,G.d);
for i=1:G.d
    max_val(i) = max(D.j(:,i))+G.dx; min_val(i) = min(D.j(:,i))-G.dx;
end
x_list = [min_val(1):G.dx:max_val(1)];
y_list = [min_val(2):G.dx:max_val(2)];
z_list = [min_val(3):G.dx:max_val(3)];
[X,Y,Z] = meshgrid(x_list, y_list, z_list);
Pfull=zeros(length(y_list),length(x_list),length(z_list));

for l=1:D.n
    x_val=D.j(l,1); y_val=D.j(l,2); z_val=D.j(l,3);
    i=find(x_list==x_val); j=find(y_list==y_val); k=find(z_list==z_val); 
    Pfull(j,i,k)=Pfull(j,i,k)+D.P(1,l);
end

figure(2*n+1)
isosurface(X,Y,Z,Pfull,0.005); 
isosurface(X,Y,Z,Pfull,0.0007); 
isosurface(X,Y,Z,Pfull,0.0001); alpha(.5),
if(n==0)
    colormap(cool);
else
    colormap(autumn);
end
drawnow;

figure(2*n+2)
isosurface(X,Y,Z,Pfull,0.0001); alpha(.5),
if(n==0)
    colormap(cool);
else
    colormap(autumn);
end
drawnow;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%