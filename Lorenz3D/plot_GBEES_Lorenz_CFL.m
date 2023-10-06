%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%%%%%% 
T=1; G.start=[-11.5; -10; 9.5]; G.unc = [0.5; 0.5; 0.5];
G.dt=.0005; G.dx=0.5; G.d=3; G.sigma=4; G.b=1; G.r=48;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%%%%%%
DATA_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Bewley/GBEES/GBEES/C++/BST/Master 3D CFL Cell-Centered/Data";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

t_G = [];
for i=0:numFiles-1
    FILE_PATH = DATA_PATH + "/pdf_0-" + num2str(i) + ".txt";
    fileID = fopen(FILE_PATH, 'r'); 
    t_G(end+1) = str2double(fgetl(fileID)); 
    fclose(fileID); 
end

t_p = linspace(0, T, 6);
indices = zeros(size(t_p));
for i = 1:numel(t_p)
    absolute_differences = abs(t_G - t_p(i));
    [~, min_index] = min(absolute_differences);
    indices(i) = min_index;
end

Y0 = G.start; tspan = [0 50]; 
options = odeset('RelTol', 1e-13); % Setting a tolerance
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,G), tspan, Y0, options);

figure(1); clf; 
plot3(Y(:,1),Y(:,2),Y(:,3),'g-','linewidth',1); view(-109,14);  hold on;
lighting phong; light('Position',[-1 0 0]); drawnow;

Y0 = G.start;
[t, Y] = ode45(@(t, Y) Lorenz3D(Y,G), t_G, Y0, options);

figure(2); clf; 
view(-109,14);  hold on;
plot3(Y(:,1),Y(:,2),Y(:,3),'k-','linewidth',2); view(-109,14);  hold on;
plot3(Y(1,1),Y(1,2),Y(1,3),'k*'),plot3(Y(end,1),Y(end,2),Y(end,3),'k*');
lighting phong; light('Position',[-1 0 0]); drawnow;

for P=1:200
    Y0=[normrnd(G.start(1),G.unc(1)); normrnd(G.start(2),G.unc(2)); normrnd(G.start(3),G.unc(3))];
    [t, Y] = ode45(@(t, Y) Lorenz3D(Y,G), t_G, Y0, options);
    plot3(Y(:,1),Y(:,2),Y(:,3),'c-.','linewidth',0.3); 
    for i=indices
        plot3(Y(i,1),Y(i,2),Y(i,3),'k+');
        drawnow;
    end
end

for i=indices
    FILE_PATH = DATA_PATH + "/pdf_0-" + num2str(i-1) + ".txt";    
    [D.P, D.j, D.n] = parseGBEES(FILE_PATH);
    Rotate_Plot(D,G);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=Lorenz3D(y,G)                          
    f=[G.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -G.b*y(3)+y(1)*y(2)-G.b*G.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rotate_Plot(D,G)        
    max_val = zeros(1,G.d); min_val = zeros(1,G.d);
    for i=1:G.d
        max_val(i) = max(D.j(:,i))+G.dx; min_val(i) = min(D.j(:,i))-G.dx;
    end
    x_list = [min_val(1):G.dx:max_val(1)];
    y_list = [min_val(2):G.dx:max_val(2)];
    z_list = [min_val(3):G.dx:max_val(3)];
    [X,Y,Z] = meshgrid(x_list, y_list, z_list);
    Pfull=zeros(length(y_list),length(x_list),length(z_list));
    
    for l=1:D.n
        x_val=D.j(l,1); y_val=D.j(l,2); z_val=D.j(l,3);
        i=find(x_list==x_val); j=find(y_list==y_val); k=find(z_list==z_val); 
        Pfull(j,i,k)=Pfull(j,i,k)+D.P(1,l);
    end
    
    figure(1)
    isosurface(X,Y,Z,Pfull,0.005); 
    isosurface(X,Y,Z,Pfull,0.0007); 
    isosurface(X,Y,Z,Pfull,0.0001); alpha(.5),
    colormap(cool);
    %set(gcf, 'Color' , 'black' );
    %set(gcf, 'NumberTitle', 'off', 'Name', 'Exploiting Sparsity'); 
    %set(gcf, 'renderer', 'Painters') 
    drawnow;
    
    figure(2)
    isosurface(X,Y,Z,Pfull,0.0001); alpha(.5),
    colormap(cool);
    %set(gcf, 'Color' , 'black' );
    %set(gcf, 'NumberTitle', 'off', 'Name', 'Exploiting Sparsity'); 
    %set(gcf, 'renderer', 'Painters') 
    drawnow;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%