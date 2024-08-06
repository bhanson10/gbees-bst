clear all; close all; clc;

% plot_PCR3BP.m
% Benjamin Hanson, 2024

%% Initial Condition
const.mu = 2.528017528540000E-5; const.LU = 668519; const.TU = 48562;  
const.U = [const.LU, const.LU, const.LU/const.TU, const.LU/const.TU]; 
const.T = 2.6513344042156235E+0; const.r =  1560.8*500;
rv.start = [1.017714765; -1.069793E-20; -1.197784E-13; 1.187104E-2]; 

%% Initializing Figures
lbl.XString = '$x$ (km)'; lbl.YString = '$y$ (km)';
initialize_figures('n', 1, 'spacing', {[50 100 700 700]},'lbl', lbl, 'axs', {'equal'}); 
set(gca, 'FontName', 'Times', 'FontSize', 14);

% Loading Europa
europa_2d = imread('./misc/Europa_2D.png');
europa_2d = imresize(europa_2d, 1.5);
[height, width, ~] = size(europa_2d);
image([(1-const.mu)*const.LU-width/2, (1-const.mu)*const.LU+width/2], [-height/2, height/2], europa_2d);

lbl.XString = '$v_x$ (km/s)';
lbl.YString = '$v_y$ (km/s)';
initialize_figures('n', 2, 'spacing', {[800 100 700 700]}, 'lbl', lbl, 'axs', {'equal'}); 
set(gca, 'FontName', 'Times', 'FontSize', 14);

%% Truth
tspan = [0,const.T]; 
x0 = rv.start;
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
[t, x] = ode87(@(t, x) PCR3BP(t, x, const), tspan, x0, options);

x(:,1:4) = x(:,1:4).*const.U; 
t = t.*const.TU; 

%% Plotting Nominal Trajectories
figure(1); 
plot(x((t <= 14*3600),1),x((t <= 14*3600),2),'r-','LineWidth',2,'DisplayName','Nominal');
plot(x(:,1),x(:,2),'r--','LineWidth',1,'DisplayName','Nominal');
drawnow;
figure(2); 
plot(x((t <= 14*3600),3),x((t <= 14*3600),4),'r-','LineWidth',2,'DisplayName','Nominal');
plot(x(:,3),x(:,4),'r--','LineWidth', 1, 'DisplayName','Nominal');
drawnow;

%% GBEES
NM = 1; 
clear p; p.name = "GBEES"; p.alpha = [0.2, 0.4, 0.6]; p.color = 'b';
P_DIR = "./c/Data/PCR3BP/PDFs/P"; % or "./c/Data/Lorenz3D/PDFs/P";

count = 1;
for nm=0:NM-1

    P_DIR_SUB = P_DIR + num2str(nm); 
    FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);
    
    for i=[0,num_files-1]
        P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";

        [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE, const);
    
        xest_gbees{count} = zeros(size(x_gbees(1,:)));
        for j=1:n_gbees
            xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
        end

        figure(1); 
        plot_nongaussian_surface(x_gbees(:,1:2),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
        figure(2); 
        plot_nongaussian_surface(x_gbees(:,3:4),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)], p);
        drawnow; 
        
        count = count + 1;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = PCR3BP(t, x, const)
    x1 = [x(3); x(4); 2*x(4)+x(1)-(const.mu*(x(1)-1+const.mu)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*(x(1)+const.mu)/(((x(1)+const.mu)^2+x(2)^2)^(1.5))); -2*x(3)+x(2)-(const.mu*x(2)/(((x(1)-1+const.mu)^2+x(2)^2)^(1.5)))-((1-const.mu)*x(2)/(((x(1)+const.mu)^2+x(2)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename, const)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        x(count, :) = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5})].*const.U';
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%