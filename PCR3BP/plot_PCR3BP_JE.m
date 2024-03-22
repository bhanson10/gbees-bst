close all; clc; clear all; 
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po = readmatrix('./Initial Conditions/periodic_orbits_JE.csv'); const.d = 4; 
const.LU = 668519; const.TU = 48562; const.mu = po(12);
rv.start=[po(2); po(3); po(5); po(6)]; n 
const.J = po(8); const.SI = po(11); 
const.T = po(9); nm = 2; size = 500;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
norm = 0; % 0: plot is in real units, 1: plot is normalized
initialize_figures(rv,const,norm); 

% Plotting MC trajectories
for k=0:nm
    DATA_PATH = append("./MC/Trajectories/Jupiter-Europa/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList); 
    
    x_mc = zeros(size,numFiles);
    y_mc = zeros(size,numFiles);
    vx_mc = zeros(size,numFiles);
    vy_mc = zeros(size,numFiles);
    for i=0:numFiles-1
        FILE_PATH = DATA_PATH + "/MC_" + num2str(i) + ".txt"; 
        fileID = fopen(FILE_PATH, 'r'); 
        t = str2double(fgetl(fileID)); 
        
        count = 1; 
        while(count <= size)
            line = split(fgetl(fileID)); % Read a line as a string
            x_mc(count, i+1) = str2double(line{1});
            y_mc(count, i+1) = str2double(line{2});
            vx_mc(count, i+1) = str2double(line{3});
            vy_mc(count, i+1) = str2double(line{4});
            count = count + 1; 
        end

        % Close the file
        fclose(fileID);
    end

    for i=1:size
        if(norm)
            figure(3);
            plot(x_mc(i,:),y_mc(i,:),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off');
            drawnow; 
            figure(4);
            plot(vx_mc(i,:),vy_mc(i,:),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off'); 
            drawnow; 
        else
            figure(3);
            plot(x_mc(i,:).*const.LU,y_mc(i,:).*const.LU,'Color', '[0.85 0.85 0.85 0.5]', 'linewidth',0.3,'HandleVisibility','off'); 
            drawnow;
            figure(4);
            plot(vx_mc(i,:).*(const.LU/const.TU),vy_mc(i,:).*(const.LU/const.TU),'Color', '[0.85 0.85 0.85 0.5]', 'linewidth',0.3,'HandleVisibility','off'); 
            drawnow; 
        end
    end
end

% Plotting nominal trajectories
Y0 = rv.start; tspan = [0 const.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) PCR3BP(Y,const), tspan, Y0, options);

if(norm)
    figure(1); 
    plot(Y(:,1),Y(:,2),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    figure(2); 
    plot(Y(:,3),Y(:,4),'k-','linewidth',1.5,'DisplayName','Nominal');
    drawnow; 
    figure(3);
    plot(Y(:,1),Y(:,2),'k-','linewidth',1.5,'DisplayName','Nominal');  %plot(Y(1,1),Y(1,2),'k*','HandleVisibility','off'),plot(Y(end,1),Y(end,2),'k*','HandleVisibility','off'); 
    drawnow; 
    figure(4); 
    plot(Y(:,3),Y(:,4),'k-','linewidth',1.5,'DisplayName','Nominal');  %plot(Y(1,3),Y(1,4),'k*','HandleVisibility','off'),plot(Y(end,3),Y(end,4),'k*','HandleVisibility','off'); 
    drawnow;
else
    figure(1); 
    plot(Y(:,1).*const.LU,Y(:,2).*const.LU,'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    figure(2);
    plot(Y(:,3).*(const.LU/const.TU),Y(:,4).*(const.LU/const.TU),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow;
    figure(3); 
    plot(Y(:,1).*const.LU,Y(:,2).*const.LU,'k-','linewidth',1.5,'DisplayName','Nominal');  %plot(Y(1,1)*const.LU,Y(1,2)*const.LU,'k*','HandleVisibility','off'),plot(Y(end,1)*const.LU,Y(end,2)*const.LU,'k*','HandleVisibility','off');
    drawnow;
    figure(4); 
    plot(Y(:,3).*(const.LU/const.TU),Y(:,4).*(const.LU/const.TU),'k-','linewidth',1.5,'DisplayName','Nominal');  %plot(Y(1,3)*(const.LU/const.TU),Y(1,4)*(const.LU/const.TU),'k*'),plot(Y(end,3)*(const.LU/const.TU),Y(end,4)*(const.LU/const.TU),'k*','HandleVisibility','off'); 
    drawnow;
end

for k=0:nm

    %Plotting GBEES PDFs
    DATA_PATH = append("./GBEES/Data/Jupiter-Europa/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList);
    
    for i=0:numFiles-1
        FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 
        fileID = fopen(FILE_PATH, 'r'); 
        t = str2double(fgetl(fileID)); 
        
        [D.P, D.j, D.n] = parseGBEES(FILE_PATH, norm, const);
        
        D.save = 1; 
        Plot_PDF(D,i,k);
        fclose(fileID);
    end
    
    % Plotting MC scatters
    DATA_PATH = append("./MC/Epochs/Jupiter-Europa/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList);
    
    for i=0:numFiles-1
        FILE_PATH = DATA_PATH + "/MC_" + num2str(i) + ".txt"; 
        fileID = fopen(FILE_PATH, 'r'); 
        t = str2double(fgetl(fileID)); 
        
        x_mc = zeros(size); y_mc = zeros(size); vx_mc = zeros(size); vy_mc = zeros(size); 
        count = 1;
        while(count <= size)
            line = split(fgetl(fileID)); % Read a line as a string
            x_mc(count) = str2double(line{1});
            y_mc(count) = str2double(line{2});
            vx_mc(count) = str2double(line{3});
            vy_mc(count) = str2double(line{4});
            count = count + 1; 
        end

        if((k==0)&&(i==0))
            if(norm)
                figure(3);
                scatter(x_mc(1),y_mc(1),1,'k','filled','DisplayName','MC'); 
                figure(4);
                scatter(vx_mc(1),vy_mc(1),1,'k','filled','DisplayName','MC'); 
            else
                figure(3);
                scatter(x_mc(1).*const.LU,y_mc(1).*const.LU,1,'k','filled','DisplayName','MC'); 
                figure(4);
                scatter(vx_mc(1).*(const.LU/const.TU),vy_mc(1).*(const.LU/const.TU),1,'k','filled','DisplayName','MC'); 
            end
        end
        if(norm)
            figure(3);
            scatter(x_mc,y_mc,1,'k','filled','HandleVisibility','off'); 
            figure(4);
            scatter(vx_mc,vy_mc,1,'k','filled','HandleVisibility','off'); 
        else
            figure(3);
            scatter(x_mc.*const.LU,y_mc.*const.LU,1,'k','filled','HandleVisibility','off'); 
            figure(4);
            scatter(vx_mc.*(const.LU/const.TU),vy_mc.*(const.LU/const.TU),1,'k','filled','HandleVisibility','off'); 
        end

        % Close the file
        fclose(fileID);
    end
end

%Adding Annotations
figure(1)
ta1 = annotation('textarrow', [0.26 0.23], [0.56 0.53]);
ta1.String = 'M_1'; 
ta1.FontName = 'Times';
ta1.FontSize = 10; 
ta1.Color = [1 0 0];  
ta2 = annotation('textarrow', [0.37 0.4], [0.72 0.76]);
ta2.String = 'M_2'; 
ta2.FontName = 'Times';
ta2.FontSize = 10; 
ta2.Color = [1 0 0];  
ta3 = annotation('textarrow', [0.36 0.39], [0.32 0.27]);
ta3.String = 'M_3'; 
ta3.FontName = 'Times';
ta3.FontSize = 10; 
ta3.Color = [1 0 0];  

figure(2)
ta4 = annotation('textarrow', [0.52 0.52], [0.75 0.8]);
ta4.String = 'M_1'; 
ta4.FontName = 'Times';
ta4.FontSize = 10; 
ta4.Color = [1 0 0];   
ta5 = annotation('textarrow', [0.59 0.63], [0.47 0.455]);
ta5.String = 'M_2'; 
ta5.FontName = 'Times';
ta5.FontSize = 10; 
ta5.Color = [1 0 0];   
ta6 = annotation('textarrow', [0.46 0.41], [0.47 0.46]);
ta6.String = 'M_3'; 
ta6.FontName = 'Times';
ta6.FontSize = 10; 
ta6.Color = [1 0 0];   

figure(3)
ta1 = annotation('textarrow', [0.26 0.23], [0.56 0.53]);
ta1.String = 'M_1'; 
ta1.FontName = 'Times';
ta1.FontSize = 10; 
ta1.Color = [1 0 0];  
ta2 = annotation('textarrow', [0.37 0.4], [0.72 0.76]);
ta2.String = 'M_2'; 
ta2.FontName = 'Times';
ta2.FontSize = 10; 
ta2.Color = [1 0 0];  
ta3 = annotation('textarrow', [0.36 0.39], [0.32 0.27]);
ta3.String = 'M_3'; 
ta3.FontName = 'Times';
ta3.FontSize = 10; 
ta3.Color = [1 0 0];  

figure(4)
ta4 = annotation('textarrow', [0.52 0.52], [0.75 0.8]);
ta4.String = 'M_1'; 
ta4.FontName = 'Times';
ta4.FontSize = 10; 
ta4.Color = [1 0 0];   
ta5 = annotation('textarrow', [0.59 0.63], [0.47 0.455]);
ta5.String = 'M_2'; 
ta5.FontName = 'Times';
ta5.FontSize = 10; 
ta5.Color = [1 0 0];   
ta6 = annotation('textarrow', [0.46 0.41], [0.47 0.46]);
ta6.String = 'M_3'; 
ta6.FontName = 'Times';
ta6.FontSize = 10; 
ta6.Color = [1 0 0];   
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

    f1 = figure(1); clf; hold all; f1.Position = [250 400 500 350]; 
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
    drawnow; 
    
    f2 = figure(2); clf; hold all; f2.Position = [750 400 500 350]; 
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
    drawnow; 

    f3 = figure(3); clf; hold all; f3.Position = [250 0 500 350];
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
    end
    drawnow; 

    f4 = figure(4); clf; hold all; f4.Position = [750 0 500 350]; 
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlim([-0.85, 0.85])
        ylim([-1.2, 1])
        xlabel("v_x (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (LU/TU)", 'FontSize', 16, 'FontName', 'Times');
    else
        xlim([-0.85*(const.LU/const.TU), 0.85*(const.LU/const.TU)])
        ylim([-1.2*(const.LU/const.TU), 1*(const.LU/const.TU)])
        xlabel("v_x (km/s)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("v_y (km/s)", 'FontSize', 16, 'FontName', 'Times');
    end
    drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_PDF(D,flag1,flag2)        
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
    
    figure(1);
    if(flag2==0) 
        if(flag1==0)
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "DisplayName", "GBEES"); 
            legend;
        else
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
        end
    elseif(flag2==1)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
    elseif(flag2==2)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
    end
    drawnow; 

    figure(3);
    if(flag2==0) 
        if(flag1==0)
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "DisplayName", "GBEES"); 
            legend;
        else
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
        end
    elseif(flag2==1)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
    elseif(flag2==2)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
    end
    drawnow; 
    
    %file_str = "pos_PDF_" + num2str(flag2) + "_" + num2str(flag1) + ".mat";
    %save(file_str, 'X', 'Y', 'pos_Pfull');

    figure(2);
    if(flag2==0) 
        if(flag1==0)
            contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "DisplayName", "GBEES"); 
            legend;
        else
            contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "HandleVisibility", "off"); 
        end
    elseif(flag2==1)
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "HandleVisibility", "off"); 
    elseif(flag2==2)
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "HandleVisibility", "off"); 
    end
    drawnow; 

    figure(4);
    if(flag2==0) 
        if(flag1==0)
            contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "DisplayName", "GBEES"); 
            legend;
        else
            contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "HandleVisibility", "off"); 
        end
    elseif(flag2==1)
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "HandleVisibility", "off"); 
    elseif(flag2==2)
        contour(VX,VY,vel_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "HandleVisibility", "off"); 
    end
    drawnow; 

    %file_str = "vel_PDF_" + num2str(flag2) + "_" + num2str(flag1) + ".mat";
    %save(file_str, 'VX', 'VY', 'vel_Pfull');
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