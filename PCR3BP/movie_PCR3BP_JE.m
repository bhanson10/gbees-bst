close all; clc; clear all; 
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po = readmatrix('./Initial Conditions/periodic_orbits_JE.csv'); const.d = 4; 
const.LU = 668519; const.TU = 48562; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; 
const.J = po(8); const.SI = po(11);  TOL = 0.0001;
const.T = po(9); nm = 2;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
norm = 0; % 0: plot is in real units, 1: plot is normalized
initialize_figures(rv,const,norm); 

% Plotting nominal trajectories
Y0 = rv.start; tspan = [0 const.T]; 
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
%for k=0:nm
for k=0:nm
    fileList = dir(fullfile("./GBEES/Movie Data/Jupiter-Europa/M" + num2str(k), '*.txt'));  % List only .txt files
    numFiles = numel(fileList);
    
    %{
    if(k==0)
        subplot(1,2,1); 
        ta1 = annotation('textarrow', [0.209259259259258,0.180555555555556], [0.543333333333337,0.508333333333333]);
        ta1.String = 'M_1'; 
        ta1.FontName = 'Times';
        ta1.FontSize = 10; 
        ta1.Color = [1 0 0];  
        subplot(1,2,2); 
        ta2 = annotation('textarrow', [0.737999999999998,0.738095238095238], [0.711666666666669,0.77]);
        ta2.String = 'M_1'; 
        ta2.FontName = 'Times';
        ta2.FontSize = 10; 
        ta2.Color = [1 0 0];   
    elseif(k==1)
        subplot(1,2,1); 
        ta1 = annotation('textarrow', [0.235449735449735,0.250661375661376], [0.671666666666667,0.73]);
        ta1.String = 'M_2'; 
        ta1.FontName = 'Times';
        ta1.FontSize = 10; 
        ta1.Color = [1 0 0];  
        subplot(1,2,2); 
        ta2 = annotation('textarrow', [0.758597883597882,0.782407407407407], [0.493333333333334,0.445]);
        ta2.String = 'M_2'; 
        ta2.FontName = 'Times';
        ta2.FontSize = 10; 
        ta2.Color = [1 0 0];   
    elseif(k==2)
        subplot(1,2,1); 
        ta1 = annotation('textarrow', [0.220899470899471,0.236772486772486], [0.305,0.261666666666667]);
        ta1.String = 'M_3'; 
        ta1.FontName = 'Times';
        ta1.FontSize = 10; 
        ta1.Color = [1 0 0];  
        subplot(1,2,2); 
        ta2 = annotation('textarrow', [0.721560846560847,0.689153439153437], [0.493333333333333,0.438333333333334]);
        ta2.String = 'M_3'; 
        ta2.FontName = 'Times';
        ta2.FontSize = 10; 
        ta2.Color = [1 0 0];   
    end
    %}

    %for i=0:numFiles-1
    for i=[0,numFiles-1]
        FILE_PATH = "./GBEES/Movie Data/Jupiter-Europa/M" + num2str(k) + "/pdf_" + num2str(i) + ".txt"; 
        fileID = fopen(FILE_PATH, 'r'); 
        t = str2double(fgetl(fileID));
        
        [D.P, D.j, D.n] = parseGBEES(FILE_PATH, norm, const);
        sgtitle(['Measurement ', num2str(k+1), ', t=', sprintf('%.*f', 4, (t+(k*2.09435))*const.TU/3600), ' hr'],'Interpreter','Latex','FontSize',16);
        Plot_PDF(D,i,k,t);
        
        frames(count) = getframe(gcf); 
        if(k==0)
            frames1(i+1) = getframe(gcf); 
        elseif(k==1)
            frames2(i+1) = getframe(gcf); 
        elseif(k==2)
            frames3(i+1) = getframe(gcf); 
        end
        subplot(1,2,1);
        delete(findobj(gca, 'Type', 'Contour'));
        subplot(1,2,2);
        delete(findobj(gca, 'Type', 'Contour'));
        fclose(fileID);
        count = count + 1; 
    end
end

%frames1(end+1) = frames2(1);
%frames2(end+1) = frames3(1); 

%create_video(frames,  'PCR3BP_JE.mp4');
%create_video(frames1, 'PCR3BP_JE_1.mp4');
%create_video(frames2, 'PCR3BP_JE_2.mp4');
%create_video(frames3, 'PCR3BP_JE_3.mp4');
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
        xlim([-1.75, 1.25])
        ylim([-1, 1])
        scatter(-const.mu,0,75,'filled','MarkerFaceColor','m','HandleVisibility','off');
        text(-const.mu,-0.08, 'Jupiter', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter(1-const.mu,0,30,'filled','MarkerFaceColor','#808080','HandleVisibility','off');
        text(1-const.mu,-0.08, 'Europa', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter(-1.00001053,0,50,'d','k','HandleVisibility','off')
        text(-0.9,-0.08, 'L_3','HorizontalAlignment','center','FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        %scatter(rv.start(1),rv.start(2),25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    else
        xlabel("x (km)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (km)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-1.75*const.LU, 1.25*const.LU])
        ylim([-1*const.LU, 1*const.LU])
        scatter(-const.mu*const.LU,0,75,'filled','MarkerFaceColor','m','HandleVisibility','off');
        text(-const.mu*const.LU,-0.08*const.LU, 'Jupiter', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter((1-const.mu)*const.LU,0,30,'filled','MarkerFaceColor','#808080','HandleVisibility','off');
        text((1-const.mu)*const.LU,-0.08*const.LU, 'Europa', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        scatter(-1.00001053*const.LU,0,50,'d','k','HandleVisibility','off')
        text(-0.9*const.LU,-0.08*const.LU, 'L_3','HorizontalAlignment','center','FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Times');
        %scatter(rv.start(1)*const.LU,rv.start(2)*const.LU,25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    end
    set(gcf, 'Color' , 'white' );
    drawnow; 
    
    subplot(1,2,2); hold on; %legend; 
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlabel("v_x (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.85, 0.85])
        ylim([-1.2, 1])
        %scatter(rv.start(3),rv.start(4),25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    else
        xlabel("v_x (km/s)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (km/s)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.85*(const.LU/const.TU), 0.85*(const.LU/const.TU)])
        ylim([-1.2*(const.LU/const.TU), 1*(const.LU/const.TU)])
        %scatter(rv.start(3)*(const.LU/const.TU),rv.start(4)*(const.LU/const.TU),25,'filled','MarkerFaceColor','k','HandleVisibility','off')
    end
    set(gcf, 'Color' , 'white' );
    drawnow; 
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
    if((abs(t-0)<1e-5)||(abs(t-0.698118)<1e-5)||(abs(t-1.39624)<1e-5)||(abs(t-2.09435)<1e-5))
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, 'HandleVisibility', 'off');
    else
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7);
    end

    drawnow; 

    subplot(1,2,2);
    if((abs(t-0)<1e-5)||(abs(t-0.698118)<1e-5)||(abs(t-1.39624)<1e-5)||(abs(t-2.09435)<1e-5))
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, 'HandleVisibility', 'off');
    else
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7);
    end
    drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, j, n] = parseGBEES(filename, norm, const)
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