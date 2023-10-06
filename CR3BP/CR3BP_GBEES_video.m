close all; clc; clear all;
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po = readmatrix('./Data/PCR3BP/periodic_orbits.csv'); 
const.LU = 389703; const.TU = 382981; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; 
const.J = po(8); const.SI = po(11); 
rv.unc = [5/const.LU; 5/const.LU; (5E-5)*const.TU/const.LU; (5E-5)*const.TU/const.LU];
%const.T = po(9); const.dt = 0.0005; const.dx = rv.unc./4; const.n = 200; 
const.T = 0.05; const.dt = 1E-5; const.dx = rv.unc./4; const.n = 10; 
const.timesteps = round(const.T/const.dt); const.FrameRate = round(5E-5/const.dt);
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
    if(mod(i,round(const.T/(100*const.dt)))==0)
        file = "./Data/PCR3BP/pdf_0-" + string(i) + ".txt";
        blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [5 inf]);
        D.P = blob(1,:)'; D.x = blob(2:3,:)'; D.v = blob(4:5,:)'; D.n = length(D.P);
        [X,Y,Z,VX,VY,VZ] = Plot_PDF(D,const,i);

        sgtitle(['iter = ', num2str(i), ', t (TU) = ', num2str(i*const.dt), ', \Delta', 't (TU) = ', num2str(const.dt)]);
        subplot(1,2,1); cla(subplot(1,2,1));
        plot(ys(1,:),ys(2,:),'g-','linewidth',1);
        scatter(ys(1,end),ys(2,end),25,'filled','MarkerFaceColor','g');
        %surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.5)
        contour(X,Y,Z,linspace(min(Z,[],'all'),max(Z,[],'all'),15));
        xlim([min(X,[],'all')-0.25*(max(X,[],'all')-min(X,[],'all')),max(X,[],'all')+0.25*(max(X,[],'all')-min(X,[],'all'))])
        ylim([min(Y,[],'all')-0.25*(max(Y,[],'all')-min(Y,[],'all')),max(Y,[],'all')+0.25*(max(Y,[],'all')-min(Y,[],'all'))])
        axis square; 

        subplot(1,2,2); cla(subplot(1,2,2)); 
        plot(ys(3,:),ys(4,:),'r-','linewidth',1);
        scatter(ys(3,end),ys(4,end),25,'filled','MarkerFaceColor','r');
        %surf(VX,VY,VZ,'EdgeColor','none','FaceAlpha',0.5)
        contour(VX,VY,VZ,linspace(min(VZ,[],'all'),max(VZ,[],'all'),15));
        xlim([min(VX,[],'all')-0.25*(max(VX,[],'all')-min(VX,[],'all')),max(VX,[],'all')+0.25*(max(VX,[],'all')-min(VX,[],'all'))])
        ylim([min(VY,[],'all')-0.25*(max(VY,[],'all')-min(VY,[],'all')),max(VY,[],'all')+0.25*(max(VY,[],'all')-min(VY,[],'all'))])
        axis square; 
        
        frames(count) = getframe(gcf); count = count + 1; 
        drawnow; 
    end
end

%create_video(frames, const, 'PCR3BP_GBEES.mp4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Z,VX,VY,VZ] = Plot_PDF(D,const,k)  
    pos = [D.x D.P]; vel = [D.v D.P]; 
    
    xyPositions = unique(pos(:, 1:2), 'rows');
    combinedList = zeros(size(xyPositions, 1), 3);
    for i = 1:size(xyPositions, 1)
        x = xyPositions(i, 1); y = xyPositions(i, 2);
        idx = find(pos(:, 1) == x & pos(:, 2) == y);
        combinedProb = sum(pos(idx, 3));
        combinedList(i, :) = [x, y, combinedProb];
    end
    pos = combinedList; 
    [X, Y] = meshgrid(min(pos(:,1)):const.dx(1):max(pos(:,1)), min(pos(:,2)):const.dx(2):max(pos(:,2)));
    Z = griddata(pos(:,1)', pos(:,2)', pos(:,3)', X, Y); 
   
    xyPositions = unique(vel(:, 1:2), 'rows');
    combinedList = zeros(size(xyPositions, 1), 3);
    for i = 1:size(xyPositions, 1)
        x = xyPositions(i, 1); y = xyPositions(i, 2);
        idx = find(vel(:, 1) == x & vel(:, 2) == y);
        combinedProb = sum(vel(idx, 3));
        combinedList(i, :) = [x, y, combinedProb];
    end
    vel = combinedList; 
    [VX, VY] = meshgrid(min(vel(:,1)):const.dx(3):max(vel(:,1)), min(vel(:,2)):const.dx(4):max(vel(:,2)));
    VZ = griddata(vel(:,1)', vel(:,2)', vel(:,3)', VX, VY); 
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
        if((i==1)||(i==const.timesteps))
            for j=1:(const.FrameRate*3)
                frame = F(i);
                writeVideo(writerObj, frame);
            end
        end
        frame = F(i);
        writeVideo(writerObj, frame);    
    end
    close(writerObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%