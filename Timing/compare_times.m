clear all; close all; clc; 

DATASET = "./Jupiter-Europa/"; TU = 48562;

f1 = fopen(DATASET + "runtime1.txt", 'r');
f2 = fopen(DATASET + "runtime2.txt", 'r');

pf = 1; % Choose which speed test to perform, 1 or 2

if(pf==1)
    [rt1, rt2, pt1, pt2, s1, s2] = parse_files_1(f1,f2);

    f1 = figure(1); clf; hold all; f1.Position = [400 250 700 475];
    l = legend; l.Location = "Northwest"; l.FontSize = 14; l.FontName = "Times"; 
    set(gca, 'FontName' , 'Times','FontSize',14);
    xlabel("Simulation time (s)", 'FontSize', 18, 'FontName', 'Times');
    xlim([rt1(1).*TU rt1(end).*TU]);
    %xlim([rt1(1) rt1(end)]);

    yyaxis left;
    plot(rt1.*TU, pt1, 'b-', 'LineWidth', 1, 'DisplayName','GBEES');
    plot(rt2.*TU, pt2, 'b--', 'LineWidth', 1,'DisplayName','MC');
    ylabel("Program time (s)", 'FontSize', 18, 'FontName', 'Times');
    
    yyaxis right;
    plot(rt1.*TU, s1, 'r-', 'LineWidth', 1, 'DisplayName','GBEES');
    plot(rt2.*TU, s2, 'r--', 'LineWidth', 1, 'DisplayName','MC');
    ylabel("Number of cells", 'FontSize', 18, 'FontName', 'Times');

elseif(pf==2)
    [rt1, rt2, pt1, pt2, s1] = parse_files_2(f1,f2);

    f1 = figure(1); clf; hold all; f1.Position = [400 250 700 475]; 
    l = legend; l.Location = "Northwest"; l.FontSize = 14; l.FontName = "Times"; 
    set(gca, 'FontName' , 'Times','FontSize',14);
    xlabel("Simulation time (TU)", 'FontSize', 18, 'FontName', 'Times');
    ylabel("Program time (s)", 'FontSize', 18, 'FontName', 'Times');
    xlim([rt1(1) rt1(end)]);

    yyaxis left;
    plot(rt1, pt1, 'b-', 'LineWidth', 1, 'DisplayName','GBEES');
    plot(rt2, pt2, 'b--', 'LineWidth', 1,'DisplayName','MC');

    yyaxis right;
    plot(rt1, s1, 'r-', 'LineWidth', 1, 'DisplayName','GBEES');
    ylabel("Number of cells", 'FontSize', 18, 'FontName', 'Times');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rt1, rt2, pt1, pt2, s1, s2] = parse_files_1(f1,f2)
    rt1 = []; rt2 = rt1; pt1 = rt1; pt2 = rt1; s1 = rt1; s2 = rt1; 
    
    while ~feof(f1)
        line = split(fgetl(f1)); 
        pt1(end+1) = str2double(line{5});
        rt1(end+1) = str2double(line{9});
        size = strsplit(line{13},'/');
        s1(end+1) = str2double(size{1}); 
    end

    while ~feof(f2)
        line = split(fgetl(f2)); 
        pt2(end+1) = str2double(line{5});
        rt2(end+1) = str2double(line{9});
        size = strsplit(line{13},'/');
        s2(end+1) = str2double(size{1}); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rt1, rt2, pt1, pt2, s1] = parse_files_2(f1,f2)
    rt1 = []; rt2 = rt1; pt1 = rt1; pt2 = rt1; s1 = rt1; 
    
    while ~feof(f1)
        line = split(fgetl(f1)); 
        pt1(end+1) = str2double(line{5});
        rt1(end+1) = str2double(line{9});
        size = strsplit(line{13},'/');
        s1(end+1) = str2double(size{1}); 
    end

    while ~feof(f2)
        line = split(fgetl(f2)); 
        pt2(end+1) = str2double(line{5});
        rt2(end+1) = str2double(line{9});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%