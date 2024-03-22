close all; clc; clear all; 
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po1 = readmatrix('./Initial Conditions/full_periodic_orbits_SE.csv'); const1.d = 4;  
const1.LU = 149597871; const1.TU = 5022635; const1.mu = po1(12);
rv1.start=[po1(2); po1(3); po1(5); po1(6)]; 
const1.T = po1(9);

po2 = readmatrix('./Initial Conditions/close-approach_periodic_orbits_SE.csv'); const.d = 4; 
const.LU = const1.LU; const.TU = const1.TU; const.mu = const1.mu; 
rv2.start=[po2(2); po2(3); po2(5); po2(6)]; 
const.T = po2(9); nm = 0; size = 500; 
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
norm = 0; % 0: plot is in real units, 1: plot is normalized
initialize_figures(rv2,const,norm); 

% Full Trajectory
Y0 = rv1.start; tspan = [0 const1.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) PCR3BP(Y,const1), tspan, Y0, options);

if(norm)
    subplot(1,2,1); 
    plot(Y(:,1),Y(:,2),'k--','linewidth',0.5,'HandleVisibility','off');
    drawnow; 
    subplot(1,2,2); 
    plot(Y(:,3),Y(:,4),'k--','linewidth',0.5,'HandleVisibility','off');
else
    subplot(1,2,1); 
    plot(Y(:,1).*const1.LU,Y(:,2).*const1.LU,'k--','linewidth',0.5,'HandleVisibility','off'); 
    drawnow; 
    subplot(1,2,2);
    plot(Y(:,3).*(const1.LU/const1.TU),Y(:,4).*(const1.LU/const1.TU),'k--','linewidth',0.5,'HandleVisibility','off'); 
    drawnow;
end

% Close Approach
Y0 = rv2.start; tspan = [0 const.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) PCR3BP(Y,const), tspan, Y0, options);

if(norm)
    subplot(1,2,1); 
    plot(Y(:,1),Y(:,2),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    subplot(1,2,2); 
    plot(Y(:,3),Y(:,4),'k-','linewidth',1.5,'DisplayName','Nominal');
    drawnow; 
else
    subplot(1,2,1); 
    plot(Y(:,1).*const.LU,Y(:,2).*const.LU,'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    subplot(1,2,2);
    plot(Y(:,3).*(const.LU/const.TU),Y(:,4).*(const.LU/const.TU),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow;
end

%Plotting GBEES PDFs
count = 1; 
for k=0:nm
    fileList = dir(fullfile("./GBEES/Movie Data/Sun-Earth/M" + num2str(k), '*.txt'));  % List only .txt files
    numFiles = numel(fileList);

    for i=0:numFiles-1
        FILE_PATH = "./GBEES/Movie Data/Sun-Earth/M" + num2str(k) + "/pdf_" + num2str(i) + ".txt"; 
        fileID = fopen(FILE_PATH, 'r'); 
        t = str2double(fgetl(fileID));
        
        [D.P, D.j, D.n] = parseGBEES(FILE_PATH, norm, const);
        sgtitle(['Measurement ', num2str(k+1), ', t=', sprintf('%.*f', 4, t*const.TU/3600), ' hr'],'Interpreter','Latex','FontSize',16);
        Plot_PDF(D,i,k,t);
        
        pause(1); 
        frames(count) = getframe(gcf); 
        subplot(1,2,1);
        delete(findobj(gca, 'Type', 'Contour'));
        subplot(1,2,2);
        delete(findobj(gca, 'Type', 'Contour'));
        fclose(fileID);
        count = count + 1; 
    end
end

create_video(frames, 'PCR3BP_SE.mp4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=PCR3BP(X,const)     
    x = X(1); y = X(2); vx = X(3); vy = X(4); 

    r1 = ((x+const.mu)^2+y^2)^(1.5);
    r2 = ((x-1+const.mu)^2+y^2)^(1.5);
    
    ax = 2*vy+x-(const.mu*(x-1+const.mu)/r2)-((1-const.mu)*(x+const.mu)/r1); 
    ay = -2*vx+y-(const.mu*y/r2)-((1-const.mu)*y/r1);

    f = [vx; vy; ax; ay]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures(rv,const,norm)

    fig = figure(1); clf; hold on; fig.Position = [0 150 1800 600]; axis tight; 
    set(gca, 'FontName' , 'Times','FontSize',12);
    sgtitle('Measurement 1, t=0.0000 hr', 'Interpreter','Latex','FontSize',16);

    subplot(1,2,1); hold on; %legend; 
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlabel("x (LU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (LU)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.06,0.02])
        ylim([-0.1,0.1])
        scatter(-const.mu,0,75,'filled','MarkerFaceColor','y','HandleVisibility','off');
        text(0,-0.02, 'Sun', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
    else
        xlabel("x (km)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (km)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.06*const.LU, 0.02*const.LU])
        ylim([-0.1*const.LU, 0.1*const.LU])
        scatter(-const.mu*const.LU,0,75,'filled','MarkerFaceColor','y','HandleVisibility','off');
        text(0,-0.02*const.LU, 'Sun', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
    end
    drawnow; 
    
    subplot(1,2,2); hold on; %legend;  
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlabel("v_x (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-4, 4])
        ylim([-6.7, -2.1])
    else
        xlabel("v_x (km/s)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (km/s)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-4*(const.LU/const.TU), 4*(const.LU/const.TU)])
        ylim([-6.7*(const.LU/const.TU), -2.1*(const.LU/const.TU)])
    end
    set(gcf, 'Color' , 'white' );
    drawnow; 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_PDF(D,flag1,flag2,t)        
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
    
    subplot(1,2,1);
    if((abs(t-0)<1e-5)||(abs(t-0.00523597)<1e-5)||(abs(t-0.0104719)<1e-5)||(abs(t-0.0157079)<1e-5)||(abs(t-0.0209439)<1e-5)||(abs(t-0.0261799)<1e-5)||(abs(t-0.0314158)<1e-5))
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, 'HandleVisibility', 'off');
    else
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7);
    end

    drawnow; 

    subplot(1,2,2);
    if((abs(t-0)<1e-5)||(abs(t-0.00523597)<1e-5)||(abs(t-0.0104719)<1e-5)||(abs(t-0.0157079)<1e-5)||(abs(t-0.0209439)<1e-5)||(abs(t-0.0261799)<1e-5)||(abs(t-0.0314158)<1e-5))
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, 'HandleVisibility', 'off');
    else
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7);
    end
    drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, j, n] = parseGBEES(filename,norm, const)
    P = []; j = [];

    fileID = fopen(filename, 'r'); line = fgetl(fileID); % Skip first line
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count) = str2double(line{1});
        if(norm)
            j(count, :) = [str2double(line{2}) str2double(line{3}) str2double(line{4}) str2double(line{5})];
        else
            j(count, :) = [str2double(line{2})*const.LU str2double(line{3})*const.LU str2double(line{4})*(const.LU/const.TU) str2double(line{5})*(const.LU/const.TU)];
        end
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_video(F, title)
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    for i=1:length(F)
        frame = F(i);
        writeVideo(writerObj, frame);    
    end
    close(writerObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%