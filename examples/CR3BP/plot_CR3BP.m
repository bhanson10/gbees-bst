clear all; close all; clc;

% plot_CR3BP.m
% Benjamin Hanson, 2024

%% Initial Condition
prop.mu = 2.528017528540000E-5; prop.LU = 668519; prop.TU = 48562; prop.sec_r =  1560.8;    
prop.U = [prop.LU, prop.LU, prop.LU, prop.LU/prop.TU, prop.LU/prop.TU, prop.LU/prop.TU]; 
prop.T = 2.6513344042156235E+0; prop.r = 1560.8;
ic = [1.0169962521469986; -1.0697953862098174E-20; -5.1360139986508944E-34; -1.3935178510137174E-14; 1.2575912545456968E-2; -3.1572204202682494E-33];

%% initializing figure
f1 = figure(1); clf; hold all; f1.Position = [50 100 1400 700];
tiledlayout(1, 2, 'TileSpacing','compact');

nexttile(1); hold all; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("$x$ (km)", "Interpreter","latex")
ylabel("$y$ (km)", "Interpreter","latex")
zlabel("$z$ (km)", "Interpreter","latex")
europa_3d = imread('Europa_map.jpg');
[X, Y, Z] = sphere(100);  
X = X.*prop.r + (1-prop.mu)*prop.LU; Y = Y.*prop.r; Z = Z.*prop.r;
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', europa_3d, 'EdgeColor', 'none', 'HandleVisibility', 'off'); alpha(0.75); 
drawnow;

nexttile(2); hold all; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("$v_x$ (km/s)", "Interpreter","latex")
ylabel("$v_y$ (km/s)", "Interpreter","latex")
zlabel("$v_z$ (km/s)", "Interpreter","latex")

%% truth
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
[t, x] = ode87(@(t, x) CR3BP(t, x, prop), [0, prop.T], ic, options);
x(:,1:6) = x(:,1:6).*prop.U; 
t = t.*prop.TU; 

%% nominal trajectory
nexttile(1); 
plot3(x((t <= 12*3600),1),x((t <= 12*3600),2),x((t <= 12*3600),3),'r-','LineWidth',2,'DisplayName','Nominal');
plot3(x(:,1),x(:,2),x(:,3),'r--','LineWidth',1,'DisplayName','Nominal');
drawnow;
nexttile(2); 
plot3(x((t <= 12*3600),4),x((t <= 12*3600),5),x((t <= 12*3600),6),'r-','LineWidth',2,'DisplayName','Nominal');
plot3(x(:,4),x(:,5),x(:,6),'r--','LineWidth', 1, 'DisplayName','Nominal');
drawnow;

%% monte carlo
n = 5000; 
P = diag([2.237500000000000E-8, 2.237500000000000E-8, 2.237500000000000E-8, 5.276700000000000E-7, 5.276700000000000E-7, 5.276700000000000E-7]);
x0_mc = mvnrnd(ic,P,n);
xf_mc = []; 
tspan = [0,12*3600/prop.TU]; 
for i=1:n
    [t, x] = ode87(@(t, x) CR3BP(t, x, prop), tspan, x0_mc(i,:), options);
    xf_mc(i,:) = x(end,:); 
end
x0_mc(:,1:6) = x0_mc(:,1:6).*prop.U;
xf_mc(:,1:6) = xf_mc(:,1:6).*prop.U;

nexttile(1);
scatter3(x0_mc(:,1), x0_mc(:,2), x0_mc(:,3), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.05);
scatter3(xf_mc(:,1), xf_mc(:,2), x0_mc(:,3), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.05);
drawnow; 

nexttile(2);
scatter3(x0_mc(:,4), x0_mc(:,5), x0_mc(:,6), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.05);
scatter3(xf_mc(:,4), xf_mc(:,5), x0_mc(:,6), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.05);
drawnow; 

%% GBEES
NM = 1; 
p.color = 'b'; p.alpha = [0.1, 0.2, 0.6]; 
P_DIR = "./results/c";

count = 1;
for nm=0:NM-1

    P_DIR_SUB = P_DIR + "/P" + num2str(nm); 
    FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);
    
    for i=[0, num_files - 1]
        P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";

        [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE, prop);
    
        xest_gbees{count} = zeros(size(x_gbees(1,:)));
        for j=1:n_gbees
            xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
        end

        nexttile(1); 
        plot_nongaussian_surface(x_gbees(:,1:3),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
        drawnow; 
        nexttile(2); 
        plot_nongaussian_surface(x_gbees(:,4:6),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)], p);
        drawnow; 
        count = count + 1;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xk = CR3BP(t, x, prop)
    xk = [x(4); x(5); x(6); 2*x(5)+x(1)-(prop.mu*(x(1)-1+prop.mu)/(((x(1)-1+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))-((1-prop.mu)*(x(1)+prop.mu)/(((x(1)+prop.mu)^2+x(2)^2+x(3)^2)^(1.5))); -2*x(4)+x(2)-(prop.mu*x(2)/(((x(1)-1+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))-((1-prop.mu)*x(2)/(((x(1)+prop.mu)^2+x(2)^2+x(3)^2)^(1.5))); -(prop.mu*x(3)/(((x(1)-1+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))-((1-prop.mu)*x(3)/(((x(1)+prop.mu)^2+x(2)^2+x(3)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename, prop)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        x(count, :) = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5});str2double(line{6});str2double(line{7})].*prop.U';
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%