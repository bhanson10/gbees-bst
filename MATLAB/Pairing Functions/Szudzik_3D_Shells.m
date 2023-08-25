clear; d = 10; a = -10;
count = 1;
for x=-d/2:d/2
    for y=-d/2:d/2
        for z=-d/2:d/2
            if (x+y+z==a)
                state = [x y z]; states{count} = state;
                if x<0, x_shift=-2*x-1; else; x_shift = 2*x; end
                if y<0, y_shift=-2*y-1; else; y_shift = 2*y; end
                if z<0, z_shift=-2*z-1; else; z_shift = 2*z; end

                if (x_shift ~= max(x_shift,y_shift))
                    xs = y_shift^2+x_shift;
                else
                    xs = x_shift^2+x_shift+y_shift;
                end

                if (xs ~= max(xs,z_shift))
                    keys(count) = z_shift^2+xs;
                else
                    keys(count) = xs^2+xs+z_shift;
                end

                count = count + 1;
            end
        end
    end
end

[keys,sortIdx] = sort(keys,'ascend');
states = states(sortIdx); num = length(states); 

x_coords = zeros(1,num); y_coords = zeros(1,num); z_coords = zeros(1,num);
for i=1:num
    state = states{i}; 
    x_coords(i) = state(1); y_coords(i) = state(2); z_coords(i) = state(3);
end

xv = linspace(min(x_coords), max(x_coords), 20);
yv = linspace(min(y_coords), max(y_coords), 20);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x_coords,y_coords,z_coords,X,Y);


figure(1); clf; grid on, hold on
scatter3(x_coords,y_coords,z_coords, 'b', 'filled', "DisplayName", "State");
plot3(x_coords,y_coords,z_coords, 'k', 'Linewidth', 1, "DisplayName", "Ascending Key Value");
surf(X,Y,Z, 'FaceAlpha',0.75, "DisplayName", "Shell"); 
title("Szudzik, x,y,z = [" + string(-d/2) + "," + string(d/2) + "], x+y+z=" +string(a), 'Interpreter', 'latex')
xticks([min(x_coords):max(x_coords)])
yticks([min(y_coords):max(y_coords)])
zticks([min(z_coords):max(z_coords)])
xlabel("x")
ylabel("y")
zlabel("z")
lgd = legend;
lgd.Location = "best";
lgd.FontSize = 10;
view(165,45)


