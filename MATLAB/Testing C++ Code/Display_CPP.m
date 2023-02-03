%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
T=1; G.thresh=0.00002; G.max=10000; G.start=[-11.5; -10; 9.5]; 
G.dt=.0005; dt=.005; G.dx=0.4; G.d=3; G.sigma=4; G.b=1; G.r=48; G.L=30; G.Y=eye(G.d,'int16'); 
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
figure(1); clf; y=G.start; ys=y; 
for timestep=1:10000
  k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
  ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
end
plot3(ys(1,:),ys(2,:),ys(3,:),'g-','linewidth',1); view(-109,14);  hold on;
lighting phong; light('Position',[-1 0 0]); drawnow;

figure(2); clf; y=G.start; ys=y; 
for timestep=1:T/dt
  k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
  ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
end
plot3(ys(1,:),ys(2,:),ys(3,:),'k-','linewidth',2); view(-109,14);  hold on;
plot3(ys(1,1),ys(2,1),ys(3,1),'k*'),plot3(ys(1,end),ys(2,end),ys(3,end),'k*');
lighting phong; light('Position',[-1 0 0]); drawnow;

for P=1:200; y=G.start+0.5*randn(3,1); ys=y;
  for timestep=1:T/dt
    k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
    ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
  end
  plot3(ys(1,:),ys(2,:),ys(3,:),'c-.','linewidth',0.3); 
  plot3(ys(1,1),ys(2,1),ys(3,1),'k+'), plot3(ys(1,end),ys(2,end),ys(3,end),'k+');
  plot3(ys(1,41),ys(2,41),ys(3,41),'k+'); plot3(ys(1,81),ys(2,81),ys(3,81),'k+');
  plot3(ys(1,121),ys(2,121),ys(3,121),'k+'); plot3(ys(1,161),ys(2,161),ys(3,161),'k+');
  drawnow;
end

for i=0:6
    file = "cpp_sim_data_v7_" + string(i) + ".txt";
    %file = "pdf_" + string(i*400) + ".txt";
    blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [4 inf]);
    D.P = blob(1,:); D.j = blob(2:4,:)'; D.n = length(D.P);
    Rotate_Plot(D,G);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(y,G)                          
f=[G.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -G.b*y(3)+y(1)*y(2)-G.b*G.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rotate_Plot(D,G)                         
N=round(2*G.L/G.dx)+1; M=(N-1)/2+1; Pfull=zeros(N,N,N);
for l=2:D.n, i=D.j(l,1)+M; j=D.j(l,2)+M; k=D.j(l,3)+M;
  if i>0 & i<=N & j>0 & j<=N & k>0 & k<=N, Pfull(j,i,k)=D.P(l); end, end
figure(1)
isosurface([-G.L:G.dx:G.L],[-G.L:G.dx:G.L],[-G.L:G.dx:G.L],Pfull,0.005); 
isosurface([-G.L:G.dx:G.L],[-G.L:G.dx:G.L],[-G.L:G.dx:G.L],Pfull,0.0007); 
isosurface([-G.L:G.dx:G.L],[-G.L:G.dx:G.L],[-G.L:G.dx:G.L],Pfull,0.0001); alpha(.5),
colormap(cool); axis([-G.L G.L -G.L G.L -G.L G.L]);
axis equal; axis([-15 15 -25 25 -30 20]); drawnow;

figure(2)
isosurface([-G.L:G.dx:G.L],[-G.L:G.dx:G.L],[-G.L:G.dx:G.L],Pfull,0.0001); alpha(.5),
colormap(cool); axis([-G.L G.L -G.L G.L -G.L G.L]);
axis equal; axis([-15 15 -25 25 -30 20]); drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%