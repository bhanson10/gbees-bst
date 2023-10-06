close all; clc; clear all; 
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po = readmatrix('./Data/periodic_orbits.csv'); const.d = 4; 
const.LU = 389703; const.TU = 382981; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; 
const.J = po(8); const.SI = po(11); 
rv.unc = [5/const.LU; 5/const.LU; (5E-5)*const.TU/const.LU; (5E-5)*const.TU/const.LU];
%const.T = po(9); const.dt = 0.001; const.n = 400; 
const.T = 3; const.dt = 1E-3; const.n = 200; 
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
    plot(ys(1,:),ys(2,:),'g-','linewidth',1); 
    drawnow; 
    figure(2); 
    plot(ys(3,:),ys(4,:),'r-','linewidth',1);
    drawnow; 
else
    figure(1); 
    plot(ys(1,:).*const.LU,ys(2,:).*const.LU,'g-','linewidth',1); 
    drawnow; 
    figure(2);
    plot(ys(3,:).*(const.LU/const.TU),ys(4,:).*(const.LU/const.TU),'r-','linewidth',1); 
    drawnow;
end

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
    plot(ys(1,:),ys(2,:),'k-','linewidth',2); plot(ys(1,1),ys(2,1),'k*'),plot(ys(1,end),ys(2,end),'k*'); 
    drawnow; 
    figure(4); 
    %xlim([-1.5E-3,1.5E-3])
    %ylim([0.756,0.7595])
    axis square;
    plot(ys(3,:),ys(4,:),'k-','linewidth',2); plot(ys(3,1),ys(4,1),'k*'),plot(ys(3,end),ys(4,end),'k*'); 
    drawnow;
else
    figure(3); 
    plot(ys(1,:).*const.LU,ys(2,:).*const.LU,'k-','linewidth',2); plot(ys(1,1)*const.LU,ys(2,1)*const.LU,'k*'),plot(ys(1,end)*const.LU,ys(2,end)*const.LU,'k*');
    drawnow;
    figure(4); 
    plot(ys(3,:).*(const.LU/const.TU),ys(4,:).*(const.LU/const.TU),'k-','linewidth',2); plot(ys(3,1)*(const.LU/const.TU),ys(4,1)*(const.LU/const.TU),'k*'),plot(ys(3,end)*(const.LU/const.TU),ys(4,end)*(const.LU/const.TU),'k*'); 
    drawnow;
end

for P=1:const.n; y=[normrnd(rv.start(1),rv.unc(1));normrnd(rv.start(2),rv.unc(2));normrnd(rv.start(3),rv.unc(3));normrnd(rv.start(4),rv.unc(4))]; ys=y;
  for timestep=1:round(const.T/const.dt)
      k1=RHS(y,const); k2=RHS(y+(const.dt/2)*k1,const); k3=RHS(y+(const.dt/2)*k2,const); k4=RHS(y+const.dt*k3,const);    
      ynew=y+(const.dt/6)*k1+(const.dt/3)*(k2+k3)+(const.dt/6)*k4; ys=[ys ynew]; y=ynew;
  end
  if(norm)
      figure(3); 
      plot(ys(1,:),ys(2,:),'c-.','linewidth',0.3); 
      plot(ys(1,1),ys(2,1),'g+'); plot(ys(1,round(end/5)),ys(2,round(end/5)),'k+'); plot(ys(1,2*round(end/5)),ys(2,2*round(end/5)),'k+');
      plot(ys(1,3*round(end/5)),ys(2,3*round(end/5)),'k+'); plot(ys(1,4*round(end/5)),ys(2,4*round(end/5)),'k+'); plot(ys(1,end),ys(2,end),'r+');
      figure(4); 
      plot(ys(3,:),ys(4,:),'c-.','linewidth',0.3); 
      plot(ys(3,1),ys(4,1),'g+'); plot(ys(3,round(end/5)),ys(4,round(end/5)),'k+'); plot(ys(3,2*round(end/5)),ys(4,2*round(end/5)),'k+');
      plot(ys(3,3*round(end/5)),ys(4,3*round(end/5)),'k+'); plot(ys(3,4*round(end/5)),ys(4,4*round(end/5)),'k+'); plot(ys(3,end),ys(4,end),'r+');
  else
      figure(3); 
      plot(ys(1,:).*const.LU,ys(2,:).*const.LU,'c-.','linewidth',0.3); 
      plot(ys(1,1)*const.LU,ys(2,1)*const.LU,'g+'); plot(ys(1,round(end/5))*const.LU,ys(2,round(end/5))*const.LU,'k+'); plot(ys(1,2*round(end/5))*const.LU,ys(2,2*round(end/5))*const.LU,'k+');
      plot(ys(1,3*round(end/5))*const.LU,ys(2,3*round(end/5))*const.LU,'k+'); plot(ys(1,4*round(end/5))*const.LU,ys(2,4*round(end/5))*const.LU,'k+'); plot(ys(1,end)*const.LU,ys(2,end)*const.LU,'r+');
      figure(4); 
      plot(ys(3,:).*(const.LU/const.TU),ys(4,:).*(const.LU/const.TU),'c-.','linewidth',0.3); 
      plot(ys(3,1)*(const.LU/const.TU),ys(4,1)*(const.LU/const.TU),'g+'); plot(ys(3,round(end/5))*(const.LU/const.TU),ys(4,round(end/5))*(const.LU/const.TU),'k+'); plot(ys(3,2*round(end/5))*(const.LU/const.TU),ys(4,2*round(end/5))*(const.LU/const.TU),'k+');
      plot(ys(3,3*round(end/5))*(const.LU/const.TU),ys(4,3*round(end/5))*(const.LU/const.TU),'k+'); plot(ys(3,4*round(end/5))*(const.LU/const.TU),ys(4,4*round(end/5))*(const.LU/const.TU),'k+'); plot(ys(3,end)*(const.LU/const.TU),ys(4,end)*(const.LU/const.TU),'r+');
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
        r1 = ((y(1,timestep)+const.mu)^2+y(2,timestep)^2)^(1/2);
        r2 = ((y(1,timestep)-1+const.mu)^2+y(2,timestep)^2)^(1/2);
        
        if(norm)
            J(P,timestep) = -2*((1/2)*(y(3,timestep)^2+y(4,timestep)^2) - (1/2)*(y(1,timestep)^2+y(2,timestep)^2)-((1-const.mu)/r1)-(const.mu/r2)); 
        else
            J(P,timestep) = (const.LU^2/const.TU^2)*(-2*((1/2)*(y(3,timestep)^2+y(4,timestep)^2) - (1/2)*(y(1,timestep)^2+y(2,timestep)^2)-((1-const.mu)/r1)-(const.mu/r2))); 
        end
    end
    if(norm)
        plot((1:round(const.T/const.dt)).*(const.dt),J(P,:));
    else
        plot((1:round(const.T/const.dt)).*(const.dt*const.TU/(86400)),J(P,:));
    end
end

for i=1:5
    file = "./Data/pdf_0-" + string(i*round(const.T/(5*const.dt))) + ".txt";
    blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [5 inf]);
    D.P = blob(1,:); D.j = blob(2:5,:)'; D.n = length(D.P);
    Plot_PDF(D);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(X,const)     
    x = X(1); y = X(2); vx = X(3); vy = X(4); 

    r1 = ((x+const.mu)^2+y^2)^(3/2);
    r2 = ((x-1+const.mu)^2+y^2)^(3/2);

    ax = 2*vy+x-((1-const.mu)*(x+const.mu)/r1)-((x-1+const.mu)*(const.mu)/r2);
    ay = -2*vx+y-((1-const.mu)*y/r1)-((const.mu)*y/r2);

    f = [vx; vy; ax; ay]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures(rv,const,norm)
    figure(1); clf; hold on; grid on; axis square; 
    if(norm)
        title(['J (LU^2/TU^2) = ', num2str(const.J), ', SI = ', num2str(const.SI), ', Period (TU) = ', num2str(const.T)])
        xlabel("x (LU)",'Interpreter','latex');
        ylabel("y (LU)",'Interpreter','latex');
        %xlim([-0.25, 1.25])
        %ylim([-0.1, 0.75])
        scatter(-const.mu,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter(1-const.mu,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
        scatter(rv.start(1),rv.start(2),25,'filled','MarkerFaceColor','g')
    else
        title(['J (km^2/s^2) = ', num2str(const.J*(const.LU^2/const.TU^2)), ', SI = ', num2str(const.SI), ', Period (days) = ', num2str(const.T*const.TU/(86400))])
        xlabel("x (km)",'Interpreter','latex');
        ylabel("y (km)",'Interpreter','latex');
        xlim([-0.25*const.LU, 1.25*const.LU])
        ylim([-0.75*const.LU, 0.75*const.LU])
        scatter(-const.mu*const.LU,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter((1-const.mu)*const.LU,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
        scatter(rv.start(1)*const.LU,rv.start(2)*const.LU,25,'filled','MarkerFaceColor','g')
    end
    
    figure(2); clf; hold on; grid on; axis square;   
    if(norm)
        title(['J (LU^2/TU^2) = ', num2str(const.J), ', SI = ', num2str(const.SI), ', Period (TU) = ', num2str(const.T)])
        xlabel("$v_x$ (LU/TU)",'Interpreter','latex');
        ylabel("$v_y$ (LU/TU)",'Interpreter','latex');
        %xlim([-1.25, 1.25])
        %ylim([-0.75, 1])
        scatter(rv.start(3),rv.start(4),25,'filled','MarkerFaceColor','r')
    else
        title(['J (km^2/s^2) = ', num2str(const.J*(const.LU^2/const.TU^2)), ', SI = ', num2str(const.SI), ', Period (days) = ', num2str(const.T*const.TU/(86400))])
        xlabel("$v_x$ (km/s)",'Interpreter','latex');
        ylabel("$v_y$ (km/s)",'Interpreter','latex');
        xlim([-1.25*(const.LU/const.TU), 1.25*(const.LU/const.TU)])
        ylim([-1.5*(const.LU/const.TU), 1*(const.LU/const.TU)])
        scatter(rv.start(3)*(const.LU/const.TU),rv.start(4)*(const.LU/const.TU),25,'filled','MarkerFaceColor','r')
    end

    figure(3); clf; hold on; grid on; axis square; 
    if(norm)
        title(['N = ', num2str(const.n), ', Init. unc. (LU) = ', num2str(rv.unc(1))])
        xlabel("x (LU)",'Interpreter','latex');
        ylabel("y (LU)",'Interpreter','latex');
        %xlim([-0.25, 1.25])
        %ylim([-0.1, 0.75])
        scatter(-const.mu,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter(1-const.mu,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
        
    else
        title(['N = ', num2str(const.n), ', Init. unc. (km) = ', num2str(rv.unc(1)*const.LU)])        
        xlabel("x (km)",'Interpreter','latex');
        ylabel("y (km)",'Interpreter','latex');
        %xlim([-0.25*const.LU, 1.25*const.LU])
        %ylim([-0.75*const.LU, 0.75*const.LU])
        scatter(-const.mu*const.LU,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
        scatter((1-const.mu)*const.LU,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
    end

    figure(4); clf; hold on; grid on; axis square;  
    if(norm)
        title(['N = ', num2str(const.n), ', Init. unc. (LU/TU) = ', num2str(rv.unc(3))])
        %xlim([-1.25, 1.25])
        %ylim([-0.75, 1])
        xlabel("$v_x$ $(LU/TU)$",'Interpreter','latex');
        ylabel("$v_y$ $(LU/TU)$",'Interpreter','latex');
    else
        title(['N = ', num2str(const.n), ', Init. unc. (km/s) = ', num2str(rv.unc(3)*(const.LU/const.TU))])
        %xlim([-1.25*(const.LU/const.TU), 1.25*(const.LU/const.TU)])
        %ylim([-1.5*(const.LU/const.TU), 1*(const.LU/const.TU)])
        xlabel("$v_x$ (km/s)",'Interpreter','latex');
        ylabel("$v_y$ (km/s)",'Interpreter','latex');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_PDF(D)        
    x_list = unique(D.j(:,1));
    y_list = unique(D.j(:,2));
    vx_list = unique(D.j(:,3));
    vy_list = unique(D.j(:,4));

    [X,Y] = meshgrid(x_list, y_list);
    pos_Pfull=zeros(length(y_list),length(x_list));
    [VX,VY] = meshgrid(vx_list, vy_list);
    vel_Pfull=zeros(length(vy_list),length(vx_list));
    
    mean_x = 0; mean_y = 0; mean_vx = 0; mean_vy = 0; weight_sum = 0; 
    for m=1:D.n
        x_val=D.j(m,1); y_val=D.j(m,2); vx_val=D.j(m,3); vy_val=D.j(m,4); 
        i=find(x_list==x_val); j=find(y_list==y_val); k=find(vx_list==vx_val); l=find(vy_list==vy_val); 
        pos_Pfull(j,i)=pos_Pfull(j,i)+D.P(1,m);
        vel_Pfull(l,k)=vel_Pfull(l,k)+D.P(1,m);

        mean_x = mean_x + x_val*D.P(1,m); 
        mean_y = mean_y + y_val*D.P(1,m); 
        mean_vx = mean_vx + vx_val*D.P(1,m); 
        mean_vy = mean_vy + vy_val*D.P(1,m); 
        weight_sum = weight_sum + D.P(1,m); 
    end
    mean_x = mean_x/weight_sum; mean_y = mean_y/weight_sum; mean_vx = mean_vx/weight_sum; mean_vy = mean_vy/weight_sum; 

    max_pos_P = max(pos_Pfull,[],'all'); pos_Pfull = pos_Pfull.*(1/max_pos_P); 
    max_vel_P = max(vel_Pfull,[],'all'); vel_Pfull = vel_Pfull.*(1/max_vel_P); 

    figure(3); 
    %xlim([0.647,0.6471])
    %ylim([-0.002,0.04])
    axis square;
    %surf(X,Y,pos_Pfull,'EdgeColor','none','FaceAlpha',0.5);
    %contour(X,Y,pos_Pfull,linspace(min(pos_Pfull,[],'all'),max(pos_Pfull,[],'all'),15));
    contour(X,Y,pos_Pfull,15);
    scatter(mean_x, mean_y, 100,'m','filled','pentagram');
    drawnow; 
    
    figure(4); 
    %xlim([-4E-4,1.5E-4])
    %ylim([0.756,0.7593])
    axis square;
    %surf(VX,VY,vel_Pfull,'EdgeColor','none','FaceAlpha',0.5);
    %contour(VX,VY,vel_Pfull,linspace(min(vel_Pfull,[],'all'),max(vel_Pfull,[],'all'),15));
    contour(VX,VY,vel_Pfull,15);
    scatter(mean_vx, mean_vy, 100,'m','filled','pentagram');
    drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%