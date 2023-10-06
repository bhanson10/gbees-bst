close all; clc; clear all; 
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po = readmatrix('./Data/periodic_orbits.csv'); 
const.LU = 389703; const.TU = 382981; const.mu = po(12);
rv.start=[po(2); po(3); po(4); po(5); po(6); po(7)]; 
const.J = po(8); const.SI = po(11); 
rv.unc = [5E-3; 5E-3; 5E-3; 1E-8; 1E-8; 1E-8];
const.T = po(9); const.dt = 0.001; const.dx = rv.unc./4; const.n = 200; G.dt = 1E-5;
%const.T = 5E-2; const.dt = 1E-5; const.dx = rv.unc./8; const.n = 200; G.dt = 1E-5;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
norm = 1; % 0: plot is in real units, 1: plot is normalized
initialize_figures(rv,const,norm); 

y=rv.start; ys=y; 
for timestep=1:(const.T/const.dt)
  k1=RHS(y,const); k2=RHS(y+(const.dt/2)*k1,const); k3=RHS(y+(const.dt/2)*k2,const); k4=RHS(y+const.dt*k3,const);    
  ynew=y+(const.dt/6)*k1+(const.dt/3)*(k2+k3)+(const.dt/6)*k4; ys=[ys ynew]; y=ynew;
end

if(norm)
    figure(1); 
    plot3(ys(1,:),ys(2,:),ys(3,:),'g-','linewidth',1); 
    drawnow; 
    figure(2); 
    plot3(ys(4,:),ys(5,:),ys(6,:),'r-','linewidth',1);
    drawnow; 
else
    figure(1); 
    plot(ys(1,:).*const.LU,ys(2,:).*const.LU,ys(3,:).*const.LU,'g-','linewidth',1); 
    drawnow; 
    figure(2);
    plot(ys(4,:).*(const.LU/const.TU),ys(5,:).*(const.LU/const.TU),ys(6,:).*(const.LU/const.TU),'r-','linewidth',1); 
    drawnow;
end
%{
y=rv.start; ys=y; 
for timestep=1:(const.T/const.dt)
  k1=RHS(y,const); k2=RHS(y+(const.dt/2)*k1,const); k3=RHS(y+(const.dt/2)*k2,const); k4=RHS(y+const.dt*k3,const);    
  ynew=y+(const.dt/6)*k1+(const.dt/3)*(k2+k3)+(const.dt/6)*k4; ys=[ys ynew]; y=ynew;
end

if(norm)
    figure(3);
    %xlim([0.6465,0.6476])
    %ylim([-0.002,0.04])
    axis square;
    plot3(ys(1,:),ys(2,:),ys(3,:),'k-','linewidth',2); plot3(ys(1,1),ys(2,1),ys(3,1),'k*'),plot3(ys(1,end),ys(2,end),ys(3,end),'k*'); 
    drawnow; 
    figure(4); 
    %xlim([-1.5E-3,1.5E-3])
    %ylim([0.756,0.7595])
    axis square;
    plot3(ys(4,:),ys(5,:),ys(6,:),'k-','linewidth',2); plot3(ys(4,1),ys(5,1),ys(6,1),'k*'),plot3(ys(4,end),ys(5,end),ys(6,end),'k*'); 
    drawnow;
else
    figure(3); 
    plot3(ys(1,:).*const.LU,ys(2,:).*const.LU,ys(3,:).*const.LU,'k-','linewidth',2); plot3(ys(1,1)*const.LU,ys(2,1)*const.LU,ys(3,1)*const.LU,'k*'),plot3(ys(1,end)*const.LU,ys(2,end)*const.LU,ys(3,end)*const.LU,'k*');
    drawnow;
    figure(4); 
    plot3(ys(4,:).*(const.LU/const.TU),ys(5,:).*(const.LU/const.TU),ys(6,:).*(const.LU/const.TU),'k-','linewidth',2); plot3(ys(4,1)*(const.LU/const.TU),ys(5,1)*(const.LU/const.TU),ys(6,1)*(const.LU/const.TU),'k*'),plot3(ys(4,end)*(const.LU/const.TU),ys(5,end)*(const.LU/const.TU),ys(6,end)*(const.LU/const.TU),'k*'); 
    drawnow;
end

for P=1:const.n; y=[normrnd(rv.start(1),rv.unc(1));normrnd(rv.start(2),rv.unc(2));normrnd(rv.start(3),rv.unc(3));normrnd(rv.start(4),rv.unc(4));normrnd(rv.start(5),rv.unc(5));normrnd(rv.start(6),rv.unc(6))]; ys=y;
  for timestep=1:round(const.T/const.dt)
      k1=RHS(y,const); k2=RHS(y+(const.dt/2)*k1,const); k3=RHS(y+(const.dt/2)*k2,const); k4=RHS(y+const.dt*k3,const);    
      ynew=y+(const.dt/6)*k1+(const.dt/3)*(k2+k3)+(const.dt/6)*k4; ys=[ys ynew]; y=ynew;
  end
  if(norm)
      figure(3); 
      plot3(ys(1,:),ys(2,:),ys(3,:),'c-.','linewidth',0.3); 
      plot3(ys(1,1),ys(2,1),ys(3,1),'g+'); plot3(ys(1,round(end/5)),ys(2,round(end/5)),ys(3,round(end/5)),'k+'); plot3(ys(1,2*round(end/5)),ys(2,2*round(end/5)),ys(3,2*round(end/5)),'k+');
      plot3(ys(1,3*round(end/5)),ys(2,3*round(end/5)),ys(3,3*round(end/5)),'k+'); plot3(ys(1,4*round(end/5)),ys(2,4*round(end/5)),ys(3,4*round(end/5)),'k+'); plot3(ys(1,end),ys(2,end),ys(3,end),'r+');
      figure(4); 
      plot3(ys(4,:),ys(5,:),ys(6,:),'c-.','linewidth',0.3); 
      plot3(ys(4,1),ys(5,1),ys(6,1),'g+'); plot3(ys(4,round(end/5)),ys(5,round(end/5)),ys(6,round(end/5)),'k+'); plot3(ys(4,2*round(end/5)),ys(5,2*round(end/5)),ys(6,2*round(end/5)),'k+');
      plot3(ys(4,3*round(end/5)),ys(5,3*round(end/5)),ys(6,3*round(end/5)),'k+'); plot3(ys(4,4*round(end/5)),ys(5,4*round(end/5)),ys(6,4*round(end/5)),'k+'); plot3(ys(4,end),ys(5,end),ys(6,end),'r+');
  else
      figure(3); 
      plot3(ys(1,:).*const.LU,ys(2,:).*const.LU,ys(3,:).*const.LU,'c-.','linewidth',0.3); 
      plot3(ys(1,1)*const.LU,ys(2,1)*const.LU,ys(3,1)*const.LU,'g+'); plot3(ys(1,round(end/5))*const.LU,ys(2,round(end/5))*const.LU,ys(3,round(end/5))*const.LU,'k+'); plot3(ys(1,2*round(end/5))*const.LU,ys(2,2*round(end/5))*const.LU,ys(3,2*round(end/5))*const.LU,'k+');
      plot3(ys(1,3*round(end/5))*const.LU,ys(2,3*round(end/5))*const.LU,ys(3,3*round(end/5))*const.LU,'k+'); plot3(ys(1,4*round(end/5))*const.LU,ys(2,4*round(end/5))*const.LU,ys(3,4*round(end/5))*const.LU,'k+'); plot3(ys(1,end).*const.LU,ys(2,end)*const.LU,ys(3,end)*const.LU,'r+');
      figure(4); 
      plot3(ys(4,:).*(const.LU/const.TU),ys(5,:).*(const.LU/const.TU),ys(6,:).*(const.LU/const.TU),'c-.','linewidth',0.3); 
      plot3(ys(4,1)*(const.LU/const.TU),ys(5,1)*(const.LU/const.TU),ys(6,1)*(const.LU/const.TU),'g+'); plot3(ys(4,round(end/5))*(const.LU/const.TU),ys(5,round(end/5))*(const.LU/const.TU),ys(6,round(end/5))*(const.LU/const.TU),'k+'); plot3(ys(4,2*round(end/5))*(const.LU/const.TU),ys(5,2*round(end/5))*(const.LU/const.TU),ys(6,2*round(end/5))*(const.LU/const.TU),'k+');
      plot3(ys(4,3*round(end/5))*(const.LU/const.TU),ys(5,3*round(end/5))*(const.LU/const.TU),ys(6,3*round(end/5))*(const.LU/const.TU),'k+'); plot3(ys(4,4*round(end/5))*(const.LU/const.TU),ys(5,4*round(end/5))*(const.LU/const.TU),ys(6,4*round(end/5))*(const.LU/const.TU),'k+'); plot3(ys(4,end)*(const.LU/const.TU),ys(5,end)*(const.LU/const.TU),ys(6,end)*(const.LU/const.TU),'r+');
  end
  orbits{P} = ys; 
end

%Jacobi Integral over time
figure(5); clf; hold on; grid on; 
title('Jacobi Constant over Time')
if(norm)
    xlabel("$t$ (TU)",'Interpreter','latex');
    xlim([1*const.dt,round(const.T/const.dt)*const.dt])
    ylabel("J (LU$^2$/TU$^2$)",'Interpreter','latex');
else
    xlabel("$t$ (days)",'Interpreter','latex');
    xlim([1*const.dt*const.TU/86400,round(const.T/const.dt)*const.dt*const.TU/86400])
    ylabel("J (km$^2$/s$^2$)",'Interpreter','latex');
end

J = zeros(const.n,round(const.T/const.dt)); 
for P=1:const.n
    y = orbits{P}; 
    for timestep=1:round(const.T/const.dt)
        r1 = ((y(1,timestep)+const.mu)^2+y(2,timestep)^2+y(3,timestep)^2)^(1/2);
        r2 = ((y(1,timestep)-1+const.mu)^2+y(2,timestep)^2+y(3,timestep)^2)^(1/2);
        
        if(norm)
            J(P,timestep) = -2*((1/2)*(y(4,timestep)^2+y(5,timestep)^2+y(6,timestep)^2) - (1/2)*(y(1,timestep)^2+y(2,timestep)^2)-((1-const.mu)/r1)-(const.mu/r2)); 
        else
            J(P,timestep) = -2*(const.LU^2/const.TU^2)*((1/2)*(y(4,timestep)^2+y(5,timestep)^2+y(6,timestep)^2) - (1/2)*(y(1,timestep)^2+y(2,timestep)^2)-((1-const.mu)/r1)-(const.mu/r2)); 
        end
    end
    if(norm)
        plot((1:round(const.T/const.dt)).*(const.dt),J(P,:));
    else
        plot((1:round(const.T/const.dt)).*(const.dt*const.TU/(86400)),J(P,:));
    end
end

for i=0:0
    file = "./Data/pdf_0-" + string(i*round(const.T/(5*G.dt))) + ".txt";
    blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [7 inf]);
    D.P = blob(1,:); D.j = blob(2:7,:)'; D.n = length(D.P);
    Plot_PDF(D,const);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(X,const)     
    x = X(1); y = X(2); z = X(3); vx = X(4); vy = X(5); vz = X(6);

    r1 = ((x+const.mu)^2+y^2+z^2)^(3/2);
    r2 = ((x-1+const.mu)^2+y^2+z^2)^(3/2);

    ax = 2*vy+x-((1-const.mu)*(x+const.mu)/r1)-((x-1+const.mu)*(const.mu)/r2);
    ay = -2*vx+y-((1-const.mu)*y/r1)-((const.mu)*y/r2);
    az = -((1-const.mu)*z/r1)-((const.mu)*z/r2);

    f = [vx; vy; vz; ax; ay; az]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures(rv,const,norm)
    figure(1); clf; hold on; grid on; axis square; view(45,30);
    if(norm)
        title(['J (LU^2/TU^2) = ', num2str(const.J), ', SI = ', num2str(const.SI), ', Period (TU) = ', num2str(const.T)])
        xlabel("x (LU)",'Interpreter','latex');
        ylabel("y (LU)",'Interpreter','latex');
        zlabel("z (LU)",'Interpreter','latex');
        %xlim([-0.25, 1.25])
        %ylim([-0.75, 0.75])
        scatter3(-const.mu,0,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter3(1-const.mu,0,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
        scatter3(rv.start(1),rv.start(2),rv.start(3),25,'filled','MarkerFaceColor','g')
    else
        title(['J (km^2/s^2) = ', num2str(const.J*(const.LU^2/const.TU^2)), ', SI = ', num2str(const.SI), ', Period (days) = ', num2str(const.T*const.TU/(86400))])
        xlabel("x (km)",'Interpreter','latex');
        ylabel("y (km)",'Interpreter','latex');
        zlabel("z (km)",'Interpreter','latex');
        %xlim([-0.25*const.LU, 1.25*const.LU])
        %ylim([-0.75*const.LU, 0.75*const.LU])
        scatter3(-const.mu*const.LU,0,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter3((1-const.mu)*const.LU,0,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
        scatter3(rv.start(1)*const.LU,rv.start(2)*const.LU,rv.start(3)*const.LU,25,'filled','MarkerFaceColor','g')
    end
    
    figure(2); clf; hold on; grid on; axis square; view(45,30);
    if(norm)
        title(['J (LU^2/TU^2) = ', num2str(const.J), ', SI = ', num2str(const.SI), ', Period (TU) = ', num2str(const.T)])
        xlabel("$v_x$ (LU/TU)",'Interpreter','latex');
        ylabel("$v_y$ (LU/TU)",'Interpreter','latex');
        zlabel("$v_z$ (LU/TU)",'Interpreter','latex');
        %xlim([-1.25, 1.25])
        %ylim([-1.5, 1])
        scatter3(rv.start(4),rv.start(5),rv.start(6),25,'filled','MarkerFaceColor','r')
    else
        title(['J (km^2/s^2) = ', num2str(const.J*(const.LU^2/const.TU^2)), ', SI = ', num2str(const.SI), ', Period (days) = ', num2str(const.T*const.TU/(86400))])
        xlabel("$v_x$ (km/s)",'Interpreter','latex');
        ylabel("$v_y$ (km/s)",'Interpreter','latex');
        zlabel("$v_z$ (km/s)",'Interpreter','latex');
        %xlim([-1.25*(const.LU/const.TU), 1.25*(const.LU/const.TU)])
        %ylim([-1.5*(const.LU/const.TU), 1*(const.LU/const.TU)])
        scatter3(rv.start(4)*(const.LU/const.TU),rv.start(5)*(const.LU/const.TU),rv.start(6)*(const.LU/const.TU),25,'filled','MarkerFaceColor','r')
    end

    figure(3); clf; hold on; grid on; axis square; view(45,30);
    if(norm)
        title(['N = ', num2str(const.n), ', Init. unc. (LU) = ', num2str(rv.unc(1))])
        xlabel("x (LU)",'Interpreter','latex');
        ylabel("y (LU)",'Interpreter','latex');
        zlabel("z (LU)",'Interpreter','latex');
        %xlim([-0.25, 1.25])
        %ylim([-0.75, 0.75])
        scatter3(-const.mu,0,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter3(1-const.mu,0,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
        
    else
        title(['N = ', num2str(const.n), ', Init. unc. (km) = ', num2str(rv.unc(1)*const.LU)])        
        xlabel("x (km)",'Interpreter','latex');
        ylabel("y (km)",'Interpreter','latex');
        zlabel("z (km)",'Interpreter','latex');
        %xlim([-0.25*const.LU, 1.25*const.LU])
        %ylim([-0.75*const.LU, 0.75*const.LU])
        scatter3(-const.mu*const.LU,0,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter3((1-const.mu)*const.LU,0,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
    end

    figure(4); clf; hold on; grid on; axis square; view(45,30);
    if(norm)
        title(['N = ', num2str(const.n), ', Init. unc. (LU/TU) = ', num2str(rv.unc(4))])
        %xlim([-1.25, 1.25])
        %ylim([-1.5, 1])
        xlabel("$v_x$ $(LU/TU)$",'Interpreter','latex');
        ylabel("$v_y$ $(LU/TU)$",'Interpreter','latex');
        zlabel("$v_z$ $(LU/TU)$",'Interpreter','latex');
    else
        title(['N = ', num2str(const.n), ', Init. unc. (km/s) = ', num2str(rv.unc(4)*(const.LU/const.TU))])
        %xlim([-1.25*(const.LU/const.TU), 1.25*(const.LU/const.TU)])
        %ylim([-1.5*(const.LU/const.TU), 1*(const.LU/const.TU)])
        xlabel("$v_x$ (km/s)",'Interpreter','latex');
        ylabel("$v_y$ (km/s)",'Interpreter','latex');
        zlabel("$v_z$ (km/s)",'Interpreter','latex');
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_PDF(D,const)        
    x_list = unique(D.j(:,1));
    y_list = unique(D.j(:,2));
    z_list = unique(D.j(:,3));
    vx_list = unique(D.j(:,4));
    vy_list = unique(D.j(:,5));
    vz_list = unique(D.j(:,6));

    [X,Y,Z] = meshgrid(x_list, y_list, z_list);
    pos_Pfull=zeros(length(y_list),length(x_list),length(z_list));
    [VX,VY,VZ] = meshgrid(vx_list, vy_list, vz_list);
    vel_Pfull=zeros(length(vy_list),length(vx_list),length(vz_list));

    for p=1:D.n
        x_val=D.j(p,1); y_val=D.j(p,2); z_val=D.j(p,3); 
        vx_val=D.j(p,4); vy_val=D.j(p,5); vz_val=D.j(p,6); 

        i=find(x_list==x_val); j=find(y_list==y_val); k=find(z_list==z_val); 
        l=find(vx_list==vx_val); m=find(vy_list==vy_val); n=find(vz_list==vz_val); 
        pos_Pfull(j,i,k)=pos_Pfull(j,i,k)+D.P(1,p);
        vel_Pfull(m,l,n)=vel_Pfull(m,l,n)+D.P(1,p);
    end

    max_pos_P = max(pos_Pfull,[],'all'); pos_Pfull = pos_Pfull.*(1/max_pos_P); 
    max_vel_P = max(vel_Pfull,[],'all'); vel_Pfull = vel_Pfull.*(1/max_vel_P); 

    figure(3); hold on; 
    xlim([x_list(1)-10*const.dx(1),x_list(end)+10*const.dx(1)])
    ylim([y_list(1)-10*const.dx(2),y_list(end)+10*const.dx(2)])
    zlim([z_list(1)-10*const.dx(3),z_list(end)+10*const.dx(3)])
    %axis equal;
    isosurface(X,Y,Z,pos_Pfull,0.01); alpha(.5),
    isosurface(X,Y,Z,pos_Pfull,0.005); alpha(.4),
    isosurface(X,Y,Z,pos_Pfull,0.0001); alpha(.3),
    drawnow; 
    
    figure(4); hold on; 
    xlim([vx_list(1)-10*const.dx(4),vx_list(end)+10*const.dx(4)])
    ylim([vy_list(1)-10*const.dx(5),vy_list(end)+10*const.dx(5)])
    zlim([vz_list(1)-10*const.dx(6),vz_list(end)+10*const.dx(6)])
    %axis equal;
    isosurface(VX,VY,VZ,vel_Pfull,0.01); alpha(.5),
    isosurface(VX,VY,VZ,vel_Pfull,0.0005); alpha(.4),
    isosurface(VX,VY,VZ,vel_Pfull,0.0001); alpha(.3),
    drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%