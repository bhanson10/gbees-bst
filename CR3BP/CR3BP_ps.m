close all; clc; clear all; 
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
po = readmatrix('./Data/periodic_orbits.csv'); 
const.LU = 389703; const.TU = 382981; const.mu = po(12);
rv.start=[po(2); po(3); po(4); po(5); po(6); po(7)]; 
const.J = po(8); const.SI = po(11); 
rv.unc = [5E-5; 5E-5; 5E-5; 5E-5; 5E-5; 5E-5];
const.T = po(9);
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%
n = 3; 
xspan = linspace(rv.start(1)-rv.unc(1), rv.start(1)+rv.unc(1),n); 
yspan = linspace(rv.start(2)-rv.unc(2), rv.start(2)+rv.unc(2),n);
zspan = linspace(rv.start(3)-rv.unc(3), rv.start(3)+rv.unc(3),n);
vxspan = linspace(rv.start(4)-rv.unc(4), rv.start(4)+rv.unc(4),n);
vyspan = linspace(rv.start(5)-rv.unc(5), rv.start(5)+rv.unc(5),n);
vzspan = linspace(rv.start(6)-rv.unc(6), rv.start(6)+rv.unc(6),n);

figure(1); clf; hold on; grid on; axis square; 
title("CR3BP Position (non-dimensionalized), Rotation-fixed Frame", 'Interpreter', 'Latex');
xlabel('x', 'Interpreter', 'Latex');
ylabel('y', 'Interpreter', 'Latex');
zlabel('z', 'Interpreter', 'Latex');
scatter3(-const.mu,0,0,50,'filled','MarkerFaceColor','b','DisplayName','Earth')
scatter3(1-const.mu,0,0,14,'filled','MarkerFaceColor','#808080','DisplayName','Moon')
view(30,30); 

figure(2); clf; hold on; grid on; axis square;  h = colorbar; 
title('Poincare Surface-of-Section $\{y=0,\dot y<0\}$','Interpreter','latex'); 
xlabel('$x$','Interpreter','latex')
ylabel('$\dot x$','Interpreter','latex')
h.Label.String = "J";

figure(3); clf; hold on; grid on; axis square; h = colorbar; 
title('Poincare Surface-of-Section $\{y=0,\dot y<0\}$','Interpreter','latex'); 
xlabel('$z$','Interpreter','latex')
ylabel('$\dot z$','Interpreter','latex')
h.Label.String = "J";

tspan = [0 20*const.T]; orbits = {}; 
for i=1:length(xspan)
    for j=1:length(yspan)
        for k=1:length(zspan)
            for l=1:length(vxspan)
                for m=1:length(vyspan)
                    for n=1:length(vzspan)
                        X0=[xspan(i); yspan(j); zspan(k); vxspan(l); vyspan(m); vzspan(n)]; 
                        r1 = ((X0(1)+const.mu)^2+X0(2)^2+X0(3)^2)^(1/2);
                        r2 = ((X0(1)-1+const.mu)^2+X0(2)^2+X0(3)^2)^(1/2);
                        orbit.J = (1/2)*(X0(4)^2+X0(5)^2+X0(6)^2) - (1/2)*(X0(1)^2+X0(2)^2)-((1-const.mu)/r1)-(const.mu/r2); 
        
                        options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @events); % Setting a tolerance
                        [t, XR, tevent, Xevent] = ode45(@(t, XR) CR3BP(XR, const), tspan, X0, options);
                        
                        orbit.X0 = X0; orbit.t = t; orbit.XR = XR; orbit.tevent = tevent; orbit.Xevent = Xevent;
                        orbits{end+1} = orbit; 
                    end
                end
            end
        end
    end
end     

for i=1:length(orbits)
    orbit = orbits{i}; 
    XR = orbit.XR; Xevent = orbit.Xevent; J = orbit.J; 
    figure(1);
    plot3(XR(:,1),XR(:,2),XR(:,3),'LineWidth', 0.5);

    if(~isempty(Xevent))
        x = Xevent(:,1); xd= Xevent(:,4); z = Xevent(:,3); zd= Xevent(:,6); 
        color = ones(size(x)).*J; 

        figure(2); 
        scatter(x, xd, 3, color, 'filled'); % Plot and use symmetry
        scatter(x,-xd, 3, color, 'filled'); % Plot and use symmetry
        figure(3); 
        scatter(z, zd, 3, color, 'filled'); % Plot and use symmetry
        scatter(z,-zd, 3, color, 'filled'); % Plot and use symmetry
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = CR3BP(X, const)
    x = X(1); y = X(2); z = X(3); vx = X(4); vy = X(5); vz = X(6);

    r1 = ((x+const.mu)^2+y^2+z^2)^(3/2);
    r2 = ((x-1+const.mu)^2+y^2+z^2)^(3/2);

    ax = 2*vy+x-((1-const.mu)*(x+const.mu)/r1)-((x-1+const.mu)*(const.mu)/r2);
    ay = -2*vx+y-((1-const.mu)*y/r1)-((const.mu)*y/r2);
    az = -((1-const.mu)*z/r1)-((const.mu)*z/r2);

    f = [vx; vy; vz; ax; ay; az]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value, isterminal, direction] = events(~, X)
    % X - Crossing Event Function
    value = X(2);

    isterminal = 0;  %   1 = end integration
                     %   0 = continue integration
    direction = -1;   %   1 = crossing with ydot > 0
                     %  -1 = crossing with ydot < 0
                     %   0 = doesn't matter (includes all)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%