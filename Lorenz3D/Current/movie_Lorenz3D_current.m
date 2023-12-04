close all; clc; clear all; 
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
const.d = 3; const.T = 1; const.dx=0.5; const.sigma=4; const.b=1; const.r=48;
rv.start=[-11.5; -10; 9.5]; rv.unc = [1; 1; 1]; nm = 0; size = 200; TOL = 0.0001;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
initialize_figures(); 

% Background Trajectory
Y0 = rv.start; tspan = [0 50]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,const), tspan, Y0, options);

subplot(1,2,1); 
plot3(Y(:,1),Y(:,2),Y(:,3),'g-','linewidth',1.5,'HandleVisibility','off');  
drawnow;

% Nominal Trajectory
Y0 = rv.start; tspan = [0 const.T]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,const), tspan, Y0, options);
subplot(1,2,1); 
plot3(Y(:,1),Y(:,2),Y(:,3),'k-','linewidth',2,'DisplayName','Nominal'); drawnow;
subplot(1,2,2); 
plot3(Y(:,1),Y(:,2),Y(:,3),'k-','linewidth',2,'DisplayName','Nominal'); drawnow;

% Plotting MC trajectories
for k=0:nm
    DATA_PATH = append("./MC/MC Freq Data/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList);
    
    x_mc = zeros(size,numFiles);
    y_mc = zeros(size,numFiles);
    z_mc = zeros(size,numFiles);
    for i=0:numFiles-1
        FILE_PATH = DATA_PATH + "/MC_" + num2str(i) + ".txt"; 
        fileID = fopen(FILE_PATH, 'r'); 
        t = str2double(fgetl(fileID)); 
        
        count = 1; 
        while (count < size)
            line = split(fgetl(fileID)); % Read a line as a string
            x_mc(count, i+1) = str2double(line{1});
            y_mc(count, i+1) = str2double(line{2});
            z_mc(count, i+1) = str2double(line{3});
            count = count + 1; 
        end

        % Close the file
        fclose(fileID);
    end

    for i=1:size
        subplot(1,2,2); 
        plot3(x_mc(i,:),y_mc(i,:),z_mc(i,:),'Color', '[0.85 0.85 0.85 0.5]','linewidth',0.3, 'HandleVisibility','off');
    end
end

dt = 1/60; 
%Plotting GBEES PDFs
for k=0:nm
    DATA_PATH = append("./GBEES/Movie Data/M", num2str(k));
    fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
    numFiles = numel(fileList);
    
    for i=0:numFiles-1
        FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 
        fileID = fopen(FILE_PATH, 'r'); 
        t = str2double(fgetl(fileID));
        
        [D.P, D.j, D.n] = parseGBEES(FILE_PATH);
        sgtitle(['Measurement 1, t=', sprintf('%.*f', 4, t), ' TU'],'Interpreter','Latex','FontSize',16);
        Plot_PDF(D,const,i,k);
        
        j = 0; 
        while(j-t <= TOL) 
            j_str = j/dt; 
            FILE_PATH = DATA_PATH + "/pdf_" + num2str(j_str) + ".txt"; 
            [D.P, D.j, D.n] = parseGBEES(FILE_PATH);
            Plot_PDF(D,const,i,k);

            if(abs(j-0) < TOL)
                subplot(1,2,1);
                text(-20,-6,8,'t = 0', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
                subplot(1,2,2);
                text(-20,-6,8,'t = 0', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
            elseif(abs(j-0.2) < TOL)
                subplot(1,2,1);
                text(-20,18,6,'t = 0.2', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
                subplot(1,2,2);
                text(-20,18,6,'t = 0.2', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
            elseif(abs(j-0.4) < TOL)
                subplot(1,2,1);
                text(-20,1,-4,'t = 0.4', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
                subplot(1,2,2);
                text(-20,1,-4,'t = 0.4', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
            elseif(abs(j-0.6) < TOL)
                subplot(1,2,1);
                text(-10,22,-1,'t = 0.6', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
                subplot(1,2,2);
                text(-10,22,-1,'t = 0.6', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
            elseif(abs(j-0.8) < TOL)
                subplot(1,2,1);
                text(-20,-8,2,'t = 0.8', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
                subplot(1,2,2);
                text(-20,-8,2,'t = 0.8', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
            elseif(abs(j-1) < TOL)
                subplot(1,2,1);
                text(10,20,15,'t = 1', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
                subplot(1,2,2);
                text(10,20,15,'t = 1', 'HorizontalAlignment','center', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Times');
            end

            j = j + 0.2; 
        end

        frames(i+1) = getframe(gcf); 
        subplot(1,2,1);
        delete(findobj(gca, 'Type', 'patch'));
        delete(findobj(gca, 'Type', 'text'));
        subplot(1,2,2);
        delete(findobj(gca, 'Type', 'patch'));
        delete(findobj(gca, 'Type', 'text'));
        fclose(fileID);
    end
end

%create_video(frames, 'test.mp4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=Lorenz3D(y,const)                          
    f=[const.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -const.b*y(3)+y(1)*y(2)-const.b*const.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialize_figures()

    fig = figure(1); clf; hold on; fig.Position = [0 150 1800 600]; axis tight; 
    set(gca, 'FontName' , 'Times','FontSize',12);
    sgtitle('Measurement 1, t=0.0000 TU', 'Interpreter','Latex','FontSize',16);

    subplot(1,2,1); hold on; 
    view(-109,14); lighting phong; light('Position',[-1 0 0]); 
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel("x", 'FontSize', 18, 'FontName', 'Times', 'Position',[-10 44 -26]);
    ylabel("y", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 -15 -42]);
    zlabel("z", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 47 8]);
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-20 20])
    xticks([-20 -10 0 10 20])
    xticklabels({'-20','-10','0','10','20'})
    ylim([-30 30])
    yticks([-30 -20 -10 0 10 20 30])
    yticklabels({'-30','-20','-10','0','10', '20', '30'})
    zlim([-30 30])
    zticks([-30 -20 -10 0 10 20 30])
    zticklabels({'-30','-20','-10','0','10', '20', '30'})
    set(gcf, 'Color' , 'white' );
    
    subplot(1,2,2); hold on; 
    view(-109,14); lighting phong; light('Position',[-1 0 0]); 
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel("x", 'FontSize', 18, 'FontName', 'Times', 'Position',[-10 44 -26]);
    ylabel("y", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 -15 -42]);
    zlabel("z", 'FontSize', 18, 'FontName', 'Times', 'Position',[0 47 8]);
    set(get(gca,'ZLabel'), 'Rotation', 0);
    xlim([-20 20])
    xticks([-20 -10 0 10 20])
    xticklabels({'-20','-10','0','10','20'})
    ylim([-30 30])
    yticks([-30 -20 -10 0 10 20 30])
    yticklabels({'-30','-20','-10','0','10', '20', '30'})
    zlim([-30 30])
    zticks([-30 -20 -10 0 10 20 30])
    zticklabels({'-30','-20','-10','0','10', '20', '30'})
    set(gcf, 'Color' , 'white' );

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_PDF(D,const,flag1,flag2)        
    max_val = zeros(1,const.d); min_val = zeros(1,const.d);
    for i=1:const.d
        max_val(i) = max(D.j(:,i))+const.dx; min_val(i) = min(D.j(:,i))-const.dx;
    end
    x_list = [min_val(1):const.dx:max_val(1)];
    y_list = [min_val(2):const.dx:max_val(2)];
    z_list = [min_val(3):const.dx:max_val(3)];
    [X,Y,Z] = meshgrid(x_list, y_list, z_list);
    Pfull=zeros(length(y_list),length(x_list),length(z_list));
    
    for l=1:D.n
        x_val=D.j(l,1); y_val=D.j(l,2); z_val=D.j(l,3);
        i=find(x_list==x_val); j=find(y_list==y_val); k=find(z_list==z_val); 
        Pfull(j,i,k)=Pfull(j,i,k)+D.P(1,l);
    end
    
    subplot(1,2,1);
    isosurface(X,Y,Z,Pfull,0.005, "HandleVisibility","off");
    isosurface(X,Y,Z,Pfull,0.0005, "HandleVisibility","off"); 
    isosurface(X,Y,Z,Pfull,0.00005, "HandleVisibility","off"); alpha(.5);
    colormap(cool);
    %set(gcf, 'NumberTitle', 'off', 'Name', 'Exploiting Sparsity'); 
    %set(gcf, 'renderer', 'Painters') 
    drawnow;
    
    subplot(1,2,2);
    isosurface(X,Y,Z,Pfull,0.00005,"HandleVisibility","off"); alpha(.25);
    colormap(cool);
    %set(gcf, 'NumberTitle', 'off', 'Name', 'Exploiting Sparsity'); 
    %set(gcf, 'renderer', 'Painters') 
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, j, n] = parseGBEES(filename)
    P = []; j = [];

    fileID = fopen(filename, 'r'); line = fgetl(fileID); % Skip first line
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count) = str2double(line{1});
        j(count, :) = [str2double(line{2}) str2double(line{3}) str2double(line{4})];
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