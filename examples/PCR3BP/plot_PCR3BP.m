clear all; close all; clc;

% plot_PCR3BP.m
% Benjamin Hanson, 2024

%% initializing system properties
prop.mu = 2.528017528540000E-5; prop.LU = 668519; prop.TU = 48562; prop.sec_r =  1560.8;  
prop.U = [prop.LU, prop.LU, prop.LU/prop.TU, prop.LU/prop.TU]; 
prop.T = 2.1988449635636802E+0;
ic = [1.0169962521469986; -1.0697953862098174E-20; -1.3935178510137174E-14; 1.2575912545456968E-2];

%% initializing figure
f1 = figure(1); clf; hold all; f1.Position = [50 100 1400 700];
tiledlayout(1, 2, 'TileSpacing','compact');

nexttile(1); hold all; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("$x$ (km)", "Interpreter","latex")
ylabel("$y$ (km)", "Interpreter","latex")
europa = nsidedpoly(1000, 'Center', [(1-prop.mu)*prop.LU, 0], 'Radius', prop.sec_r);
plot(europa, 'FaceColor', 'm', 'EdgeColor','none');

nexttile(2); hold all; 
axis("equal")
set(gca, 'FontName', 'Times', 'FontSize', 14);
xlabel("$v_x$ (km/s)", "Interpreter","latex")
ylabel("$v_y$ (km/s)", "Interpreter","latex")
 
%% truth
options = odeset('MaxStep', 1E-3, 'InitialStep', 1E-3, 'RelTol', 1e-6);
[t, x] = ode87(@(t, x) PCR3BP(t, x, prop), [0,prop.T], ic, options);
x(:,1:4) = x(:,1:4).*prop.U; 
t = t.*prop.TU; 

%% nominal trajectory
nexttile(1); 
plot(x((t <= 12*3600),1),x((t <= 12*3600),2),'r-','LineWidth',2,'HandleVisibility','off');
plot(x(:,1),x(:,2),'r--','LineWidth',1,'HandleVisibility','off');
drawnow;
nexttile(2); 
plot(x((t <= 12*3600),3),x((t <= 12*3600),4),'r-','LineWidth',2,'HandleVisibility','off');
plot(x(:,3),x(:,4),'r--','LineWidth', 1,'HandleVisibility','off');
drawnow;

%% monte carlo
n = 5000; 
P = diag([2.237500000000000E-8, 2.237500000000000E-8, 5.276700000000000E-7, 5.276700000000000E-7]);
x0_mc = mvnrnd(ic,P,n);
xf_mc = []; 
tspan = [0,12*3600/prop.TU]; 
for i=1:n
    [t, x] = ode87(@(t, x) PCR3BP(t, x, prop), tspan, x0_mc(i,:), options);
    xf_mc(i,:) = x(end,:); 
end
x0_mc(:,1:4) = x0_mc(:,1:4).*prop.U;
xf_mc(:,1:4) = xf_mc(:,1:4).*prop.U;

nexttile(1);
scatter(x0_mc(:,1), x0_mc(:,2), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
scatter(xf_mc(:,1), xf_mc(:,2), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
drawnow; 

nexttile(2);
scatter(x0_mc(:,3), x0_mc(:,4), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
scatter(xf_mc(:,3), xf_mc(:,4), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
drawnow; 

%% GBEES
NM = 1; 
p.color = 'b'; p.alpha = [0.2, 0.4, 0.6]; 
P_DIR = "<path_to_pdf>";

count = 1;
for nm=0:NM-1

    P_DIR_SUB = P_DIR + "/P" + num2str(nm); 
    FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
    num_files = numel(FILE_LIST);
    
    for i=[0,num_files-1]
        P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";

        [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE, prop);
    
        xest_gbees{count} = zeros(size(x_gbees(1,:)));
        for j=1:n_gbees
            xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
        end

        nexttile(1); 
        plot_nongaussian_surface(x_gbees(:,1:2),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
        nexttile(2); 
        plot_nongaussian_surface(x_gbees(:,3:4),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)], p);
        drawnow; 
        
        count = count + 1;
        p.display = 0; 
    end
end
nexttile(1);
ylim([-8000,8000]);

% %% GBEES - Monte
% NM = 1; 
% clear p; p.color = 'g'; p.alpha = [0.2, 0.4, 0.6];
% P_DIR = "./results/monte/P"; 
% 
% count = 1;
% for nm=0:NM-1
% 
%     P_DIR_SUB = P_DIR + num2str(nm); 
%     FILE_LIST = dir(fullfile(P_DIR_SUB, '*.txt'));  % List only .txt files
%     num_files = numel(FILE_LIST);
% 
%     for i=[0,num_files-1]
%         P_FILE = P_DIR_SUB + "/pdf_" + num2str(i) + ".txt";
% 
%         [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(P_FILE, prop);
% 
%         xest_gbees{count} = zeros(size(x_gbees(1,:)));
%         for j=1:n_gbees
%             xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
%         end
% 
%         nexttile(1); 
%         plot_nongaussian_surface(x_gbees(:,1:2),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)],p);
%         nexttile(2); 
%         plot_nongaussian_surface(x_gbees(:,3:4),P_gbees,[normpdf(1)/normpdf(0), normpdf(2)/normpdf(0), normpdf(3)/normpdf(0)], p);
%         drawnow; 
% 
%         count = count + 1;
%         p.display = 0;
%     end
% end
% 
% nexttile(1);
% ylim([-8000,8000]);
% 
% clear L; clear LH; 
% LH(1) = fill(nan, nan, nan, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'none');
% L{1} = "{ }Analytical{           }";
% LH(2) = fill(nan, nan, nan, 'FaceAlpha', 0.5, 'FaceColor', 'g', 'EdgeColor', 'none');
% L{2} = "{ }Monte{  }";
% leg = legend(LH, L, 'Orientation', 'Horizontal', 'FontSize', 18, 'FontName', 'times');
% leg.Layout.Tile = 'south';
% drawnow; 
% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = PCR3BP(t, x, prop)
    x1 = [x(3); x(4); 2*x(4)+x(1)-(prop.mu*(x(1)-1+prop.mu)/(((x(1)-1+prop.mu)^2+x(2)^2)^(1.5)))-((1-prop.mu)*(x(1)+prop.mu)/(((x(1)+prop.mu)^2+x(2)^2)^(1.5))); -2*x(3)+x(2)-(prop.mu*x(2)/(((x(1)-1+prop.mu)^2+x(2)^2)^(1.5)))-((1-prop.mu)*x(2)/(((x(1)+prop.mu)^2+x(2)^2)^(1.5)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename, prop)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        x(count, :) = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5})].*prop.U';
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%