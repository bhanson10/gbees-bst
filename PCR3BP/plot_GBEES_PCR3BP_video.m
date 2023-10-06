close all; clc; clear all;
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po = readmatrix('./Data/periodic_orbits.csv'); 
const.LU = 389703; const.TU = 382981; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; 
const.J = po(8); const.SI = po(11); 
rv.unc = [5/const.LU; 5/const.LU; (5E-5)*const.TU/const.LU; (5E-5)*const.TU/const.LU];
%const.T = po(9); const.dt = 0.0005;
const.T = 3; const.dt = 1E-3; 
const.timesteps = round(const.T/const.dt); const.FrameRate = round(5E-3/const.dt);
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
fig = figure(1); clf; hold on; fig.Position = [150 150 1200 600];
subplot(1,2,1); hold on; grid on; 
xlabel('x','FontSize',16,'Interpreter','latex');
ylabel('y','FontSize',16,'Interpreter','latex');
colorbar; 
subplot(1,2,2);  hold on; grid on; 
xlabel('$v_x$','FontSize',16,'Interpreter','latex');
ylabel('$v_y$','FontSize',16,'Interpreter','latex');
colorbar; 

y=rv.start; ys=y; count = 1; 
for i=0:const.timesteps
    if(i~=0)
        k1=RHS(y,const); k2=RHS(y+(const.dt/2)*k1,const); k3=RHS(y+(const.dt/2)*k2,const); k4=RHS(y+const.dt*k3,const);    
        ynew=y+(const.dt/6)*k1+(const.dt/3)*(k2+k3)+(const.dt/6)*k4; ys=[ys ynew]; y=ynew;
    end
    if((mod(i,round(const.T/(100*const.dt)))==0))
        file = "./Data/pdf_0-" + string(i) + ".txt";
        blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [5 inf]);
        D.P = blob(1,:); D.j = blob(2:5,:)'; D.n = length(D.P);
        [X,Y,pos_Pfull,VX,VY,vel_Pfull,mean] = Plot_PDF(D);

        sgtitle(['iter = ', num2str(i), ', t (TU) = ', num2str(i*const.dt), ', \Delta', 't (TU) = ', num2str(const.dt)]);
        subplot(1,2,1); cla(subplot(1,2,1));
        plot(ys(1,:),ys(2,:),'g-','linewidth',1);
        scatter(ys(1,end),ys(2,end),25,'filled','MarkerFaceColor','g');
        %surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.5)
        %contour(X,Y,pos_Pfull,linspace(min(pos_Pfull,[],'all'),max(pos_Pfull,[],'all'),15));
        contour(X,Y,pos_Pfull,15);
        scatter(mean.x, mean.y, 100,'m','filled','pentagram');
        xlim([min(X,[],'all')-0.25*(max(X,[],'all')-min(X,[],'all')),max(X,[],'all')+0.25*(max(X,[],'all')-min(X,[],'all'))])
        ylim([min(Y,[],'all')-0.25*(max(Y,[],'all')-min(Y,[],'all')),max(Y,[],'all')+0.25*(max(Y,[],'all')-min(Y,[],'all'))])
        axis square; 

        subplot(1,2,2); cla(subplot(1,2,2)); 
        plot(ys(3,:),ys(4,:),'r-','linewidth',1);
        scatter(ys(3,end),ys(4,end),25,'filled','MarkerFaceColor','r');
        %surf(VX,VY,VZ,'EdgeColor','none','FaceAlpha',0.5)
        %contour(VX,VY,vel_Pfull,linspace(min(vel_Pfull,[],'all'),max(vel_Pfull,[],'all'),15));
        contour(VX,VY,vel_Pfull,15);
        scatter(mean.vx, mean.vy, 100,'m','filled','pentagram');
        xlim([min(VX,[],'all')-0.25*(max(VX,[],'all')-min(VX,[],'all')),max(VX,[],'all')+0.25*(max(VX,[],'all')-min(VX,[],'all'))])
        ylim([min(VY,[],'all')-0.25*(max(VY,[],'all')-min(VY,[],'all')),max(VY,[],'all')+0.25*(max(VY,[],'all')-min(VY,[],'all'))])
        axis square; 
        
        set(gcf, 'Color', 'w');
        frames(count) = getframe(gcf); count = count + 1; 
        drawnow; 
    end
end

%create_video(frames, const, 'PCR3BP_GBEES.mp4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,pos_Pfull,VX,VY,vel_Pfull,mean] = Plot_PDF(D)  
    x_list = unique(D.j(:,1));
    y_list = unique(D.j(:,2));
    vx_list = unique(D.j(:,3));
    vy_list = unique(D.j(:,4));

    [X,Y] = meshgrid(x_list, y_list);
    pos_Pfull=zeros(length(y_list),length(x_list));
    [VX,VY] = meshgrid(vx_list, vy_list);
    vel_Pfull=zeros(length(vy_list),length(vx_list));
    
    mean.x = 0; mean.y = 0; mean.vx = 0; mean.vy = 0; weight_sum = 0; 
    for m=1:D.n
        x_val=D.j(m,1); y_val=D.j(m,2); vx_val=D.j(m,3); vy_val=D.j(m,4); 
        i=find(x_list==x_val); j=find(y_list==y_val); k=find(vx_list==vx_val); l=find(vy_list==vy_val); 
        pos_Pfull(j,i)=pos_Pfull(j,i)+D.P(1,m);
        vel_Pfull(l,k)=vel_Pfull(l,k)+D.P(1,m);

        mean.x = mean.x + x_val*D.P(1,m); 
        mean.y = mean.y + y_val*D.P(1,m); 
        mean.vx = mean.vx + vx_val*D.P(1,m); 
        mean.vy = mean.vy + vy_val*D.P(1,m); 
        weight_sum = weight_sum + D.P(1,m); 
    end
    mean.x = mean.x/weight_sum; mean.y = mean.y/weight_sum; mean.vx = mean.vx/weight_sum; mean.vy = mean.vy/weight_sum; 
end
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
function create_video(F, const, title)
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = const.FrameRate;
    open(writerObj);
    for i=1:length(F)
        frame = F(i);
        writeVideo(writerObj, frame);    
    end
    close(writerObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%