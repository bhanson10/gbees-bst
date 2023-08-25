clear; l = 0; u = 1;
count = 1;
for x=l:u
    for y=l:u
        for z=l:u
            state = [x y z]; states{count} = state;
            if x<0, x_shift=-2*x-1; else; x_shift = 2*x; end
            if y<0, y_shift=-2*y-1; else; y_shift = 2*y; end
            if z<0, z_shift=-2*z-1; else; z_shift = 2*z; end
            xs = (1/2)*(x_shift+y_shift)*(x_shift+y_shift+1)+y_shift;
            keys(count) = (1/2)*(xs+z_shift)*(xs+z_shift+1)+z_shift;
            count = count + 1;
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

% Create file name variable
filename = 'cantor-3d-10.gif';

figure; clf; grid on, hold on
s = scatter3(x_coords(1),y_coords(1),z_coords(1), 'b', 'filled', "DisplayName", "State");
p = plot3(x_coords(1),y_coords(1),z_coords(1), 'k', 'Linewidth', 1, "DisplayName", "Ascending Key Value");
title("Cantor Ascending Key, x,y,z = [" +string(l) + "," + string(u) + "]", 'Interpreter', 'latex')
xlabel("x")
ylabel("y")
zlabel("z")
xlim([l u])
ylim([l u])
zlim([l u])
view(165,45)

pause(3);
for k = 2:num

    % Updating the line
    p.XData = x_coords(1:k);
    p.YData = y_coords(1:k);
    p.ZData = z_coords(1:k);
    % Updating the point
    s.XData = x_coords(k); 
    s.YData = y_coords(k);
    s.ZData = z_coords(k);
    

    % Updating the title
    title("[x,y,z] = [" + x_coords(k) + "," + y_coords(k) + "," + z_coords(k) +"], Key = " + keys(k), 'Interpreter', 'latex')

    % Delay
    pause(0.01)
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 2
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',(0.5));
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',(0.5));
    end
end