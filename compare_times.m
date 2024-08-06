clear all; close all; clc; 

mkrs = {"-o", "-square", "->", "-pentagram"}; 
DIR = {"c","python"};
MDL = "Lorenz3D";

TU = 1; n = 4; 
names = {'c implementation',...
         'python extension'};
time_unit = '(TU)';

f = {}; 
for i=1:length(names)
    f{end+1} = fopen("./" + DIR{i} + "/Data/" + MDL + "/runtime.txt", 'r');
end

[rt, pt, s] = parse_time_files(f);

f1 = figure(1); clf; f1.Position = [200 350 700 475]; ax = axes; 
l = legend; l.Location = "Northwest"; l.FontSize = 14; l.FontName = "Times"; 
set(ax, 'FontName' , 'Times','FontSize',14);
xlabel("Simulation time " + time_unit, 'FontSize', 18, 'FontName', 'Times');

yyaxis left; hold on; 
ax.YAxis(1).Color = [0 0 1];
for i=1:length(f)
    rti = rt{i}; pti = pt{i}; 
    plot(rti.*TU, pti, mkrs{i}, 'Color', 'blue', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', names{i});
end
ylabel("Program time (s)", 'FontSize', 18, 'FontName', 'Times');

yyaxis right; hold on; 
ax.YAxis(2).Color = [1 0 0];
for i=1:length(f)
    rti = rt{i}; si = s{i}; 
    plot(rti.*TU, si, mkrs{i}, 'Color', 'red', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', names{i});
end
ylabel("Number of cells", 'FontSize', 18, 'FontName', 'Times');

% Normalized time
f2 = figure(2); clf; f2.Position = [900 350 700 475]; ax = axes; 
l = legend; l.Location = "Northeast"; l.FontSize = 14; l.FontName = "Times"; 
set(ax, 'FontName' , 'Times','FontSize',14);
xlabel("Simulation time " + time_unit, 'FontSize', 18, 'FontName', 'Times');

yyaxis left; hold on; 
rt1 = rt{1}; pt1 = pt{1}; name1 = names{1}; 
for i=2:length(f)
    pti = pt{i}; name = names{i}; 
    pti_norm = pti./pt1; 
    plot(rt1(2:end).*TU, pti_norm(2:end), mkrs{i}, 'Color', 'blue', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', append(name, ' / ', name1));
end
ylabel("Normalized program time", 'FontSize', 18, 'FontName', 'Times');

yyaxis right; hold on; 
s1 = s{1}; 
for i=2:length(f)
    si = s{i}; name = names{i}; 
    si_norm = si./s1; 
    plot(rt1(2:end).*TU, si_norm(2:end), mkrs{i}, 'Color', 'red', 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName', append(name, ' / ', name1));
end
ylabel("Normalized cell number", 'FontSize', 18, 'FontName', 'Times');

ax.YAxis(1).Color = [0 0 1];
ax.YAxis(2).Color = [1 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rt, pt, s] = parse_time_files(f)
    rt = {}; pt = {}; s = {}; 
    for i=1:length(f)
        rti = []; pti = rti; si = rti; fi = f{i};
        
        while ~feof(fi)
            line = split(fgetl(fi)); 
            pti(end+1) = str2double(line{5});
            rti(end+1) = str2double(line{9});
            size = strsplit(line{13},'/');
            si(end+1) = str2double(size{1}); 
        end
        
        rt{end+1} = rti;
        pt{end+1} = pti;
        s{end+1}  = si;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%