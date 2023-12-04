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

% Plotting MC trajectories
for k=0:nm
    DATA_PATH = append("./MC/Movie Data/Sun-Earth/M", num2str(k));
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
            plot(x_mc(i,:),y_mc(i,:),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off'); 
        else
            figure(3);
            plot(x_mc(i,:).*const.LU,y_mc(i,:).*const.LU,'Color', '[0.85 0.85 0.85 0.5]', 'linewidth',0.3,'HandleVisibility','off'); 
            drawnow;
            figure(4);
            plot(x_mc(i,:).*const.LU,y_mc(i,:).*const.LU,'Color', '[0.85 0.85 0.85 0.5]', 'linewidth',0.3,'HandleVisibility','off'); 
        end
    end
end
drawnow; 

% Full Trajectory
Y0 = rv1.start; tspan = [0 const1.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) PCR3BP(Y,const1), tspan, Y0, options);

if(norm) 
    plot(Y(:,1),Y(:,2),'k--','linewidth',0.5,'HandleVisibility','off');
    drawnow; 
    figure(4); 
    plot(Y(:,1),Y(:,2),'k--','linewidth',0.5,'HandleVisibility','off');
    drawnow; 
else
    figure(3); 
    plot(Y(:,1).*const1.LU,Y(:,2).*const1.LU,'k--','linewidth',0.5,'HandleVisibility','off'); 
    drawnow; 
    figure(4);
    plot(Y(:,1).*const1.LU,Y(:,2).*const1.LU,'k--','linewidth',0.5,'HandleVisibility','off'); 
    drawnow;
end

% Close Approach
Y0 = rv2.start; tspan = [0 const.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) PCR3BP(Y,const), tspan, Y0, options);

if(norm)
    figure(3); 
    plot(Y(:,1),Y(:,2),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    figure(4); 
    plot(Y(:,1),Y(:,2),'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
else
    figure(3); 
    plot(Y(:,1).*const.LU,Y(:,2).*const.LU,'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow; 
    figure(4);
    plot(Y(:,1).*const.LU,Y(:,2).*const.LU,'k-','linewidth',1.5,'DisplayName','Nominal'); 
    drawnow;
end

for k=0:nm
    %Plotting GBEES PDFs
    DATA_PATH = append("./GBEES/Data/Sun-Earth/ND/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList);

    FILE_PATH = DATA_PATH + "/pdf_" + num2str(numFiles-1) + ".txt"; 
    fileID = fopen(FILE_PATH, 'r'); 
    t = str2double(fgetl(fileID)); 
    
    [D1.P, D1.j, D1.n] = parseGBEES(FILE_PATH, norm, const);

    DATA_PATH = append("./GBEES/Data/Sun-Earth/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList);
    

    FILE_PATH = DATA_PATH + "/pdf_" + num2str(numFiles-1) + ".txt"; 
    fileID = fopen(FILE_PATH, 'r'); 
    t = str2double(fgetl(fileID)); 
    
    [D2.P, D2.j, D2.n] = parseGBEES(FILE_PATH, norm, const);

    Plot_PDF(D1,D2,0,k);
    fclose(fileID);

    %Plotting MC dots
    DATA_PATH = append("./MC/Data/Sun-Earth/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList);
    
    FILE_PATH = DATA_PATH + "/MC_" + num2str(numFiles-1) + ".txt"; 
    fileID = fopen(FILE_PATH, 'r'); 
    t = str2double(fgetl(fileID)); 
    
    x_mc = zeros(size); y_mc = zeros(size); 
    count = 1;
    while(count <= size)
        line = split(fgetl(fileID)); % Read a line as a string
        x_mc(count) = str2double(line{1});
        y_mc(count) = str2double(line{2});
        count = count + 1; 
    end

    if(norm)
        figure(3);
        scatter(x_mc(1),y_mc(1),2,'k','filled','DisplayName','MC');
        figure(4);
        scatter(x_mc(1),y_mc(1),2,'k','filled','DisplayName','MC');
    else
        figure(3);
        scatter(x_mc(1).*const.LU,y_mc(1).*const.LU,2,'k','filled','DisplayName','MC');
        figure(4);
        scatter(x_mc(1).*const.LU,y_mc(1).*const.LU,2,'k','filled','DisplayName','MC');
    end

    if(norm)
        figure(3);
        scatter(x_mc,y_mc,2,'k','filled','HandleVisibility','off');
        figure(4);
        scatter(x_mc,y_mc,2,'k','filled','HandleVisibility','off');
    else
        figure(3);
        scatter(x_mc.*const.LU,y_mc.*const.LU,2,'k','filled','HandleVisibility','off');
        figure(4);
        scatter(x_mc.*const.LU,y_mc.*const.LU,2,'k','filled','HandleVisibility','off');
    end

    % Close the file
    fclose(fileID);
end
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

    f3 = figure(3); clf; hold all; f3.Position = [250 200 500 350]; legend;
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlabel("x (LU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (LU)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.0241,-0.0106])
        ylim([-0.0835,-0.0748])
    else
        xlabel("x (km)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (km)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.0241*const.LU, -0.0106*const.LU])
        ylim([-0.0835*const.LU, -0.0748*const.LU])
    end
    drawnow; 

    f4 = figure(4); clf; hold all; f4.Position = [750 200 500 350]; legend;
    set(gca, 'FontName' , 'Times','FontSize',12);
    if(norm)
        xlabel("x (LU)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (LU)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.0241,-0.0106])
        ylim([-0.0835,-0.0748])
    else
        xlabel("x (km)", 'FontSize', 16, 'FontName', 'Times');
        ylabel("y (km)", 'FontSize', 16, 'FontName', 'Times');
        xlim([-0.0241*const.LU, -0.0106*const.LU])
        ylim([-0.0835*const.LU, -0.0748*const.LU])
    end
    drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_PDF(D1,D2,flag1,flag2)        
    x_list = unique(D1.j(:,1));
    y_list = unique(D1.j(:,2));
    
    disp("No diffusion:");
    disp("x:")
    disp(x_list(1))
    disp(x_list(end))
    disp("y:")
    disp(y_list(1))
    disp(y_list(end))

    [X,Y] = meshgrid(x_list, y_list);
    pos_Pfull=zeros(length(y_list),length(x_list));
    
    mean_x = 0; mean_y = 0; mean_vx = 0; mean_vy = 0; weight_sum = 0; 
    for m=1:D1.n
        x_val=D1.j(m,1); y_val=D1.j(m,2);
        i=find(x_list==x_val); j=find(y_list==y_val);
        pos_Pfull(j,i)=pos_Pfull(j,i)+D1.P(1,m);

        mean_x = mean_x + x_val*D1.P(1,m); 
        mean_y = mean_y + y_val*D1.P(1,m); 
        weight_sum = weight_sum + D1.P(1,m); 
    end
    mean_x = mean_x/weight_sum; mean_y = mean_y/weight_sum;

    max_pos_P = max(pos_Pfull,[],'all'); pos_Pfull = pos_Pfull.*(1/max_pos_P); 

    figure(3); 
    if(flag2==0) 
        if(flag1==0)
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "DisplayName", "GBEES, Q=0"); 
        else
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
        end
    elseif(flag2==1)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
    elseif(flag2==2)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
    end
    drawnow; 

    x_list = unique(D2.j(:,1));
    y_list = unique(D2.j(:,2));
    
    disp("Diffusion:");
    disp("x:")
    disp(x_list(1))
    disp(x_list(end))
    disp("y:")
    disp(y_list(1))
    disp(y_list(end))

    [X,Y] = meshgrid(x_list, y_list);
    pos_Pfull=zeros(length(y_list),length(x_list));
    
    mean_x = 0; mean_y = 0; mean_vx = 0; mean_vy = 0; weight_sum = 0; 
    for m=1:D2.n
        x_val=D2.j(m,1); y_val=D2.j(m,2);
        i=find(x_list==x_val); j=find(y_list==y_val);
        pos_Pfull(j,i)=pos_Pfull(j,i)+D2.P(1,m);

        mean_x = mean_x + x_val*D2.P(1,m); 
        mean_y = mean_y + y_val*D2.P(1,m); 
        weight_sum = weight_sum + D2.P(1,m); 
    end
    mean_x = mean_x/weight_sum; mean_y = mean_y/weight_sum;

    max_pos_P = max(pos_Pfull,[],'all'); pos_Pfull = pos_Pfull.*(1/max_pos_P); 

    figure(4); 
    if(flag2==0) 
        if(flag1==0)
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on', 'FaceAlpha', 0.7, "DisplayName", "GBEES, Q≠0"); 
        else
            contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
        end
    elseif(flag2==1)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
    elseif(flag2==2)
        contour(X,Y,pos_Pfull,6, 'LineWidth',1, 'Fill', 'on',  'FaceAlpha', 0.7, "HandleVisibility","off");
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