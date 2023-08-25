%Constant State
dimension = [2:6];
CantorCoords1D = zeros(1,5);
SzudzikCoords1D = zeros(1,5);
for i=1:5
    state = ones(1,i+1); shift_state = ShiftState(state,i+1); 
    CantorCoords1D(i) = CantorPair(shift_state);
    SzudzikCoords1D(i) = SzudzikPair(shift_state);
end

figure(1); clf; h=gca; set(h,'yscale','log'); hold on
plot(dimension,CantorCoords1D, 'r', 'Linewidth', 1, 'DisplayName', 'Cantor Pairing');
plot(dimension,SzudzikCoords1D, 'b', 'Linewidth', 1, 'DisplayName', 'Szudzik Pairing');
title('Mag. of Key vs. Dim. of State', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "northwest";
lgd.FontSize = 10;
xticks(dimension)
xlabel('Dimension of State', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('Magnitude of Key', 'Interpreter', 'Latex', 'FontSize', 10)

%Manipulating the Magnitude of the State
CantorCoords=zeros(500,5);
SzudzikCoords=zeros(500,5);

state_mag = [];
for i=2:6
    for z=1:50:5000
        state = z*ones(1,i); shift_state = ShiftState(state,i); 
        state_norm = norm(shift_state); state_mag = [state_mag state_norm];
    end
end

[X,Y] = meshgrid(2:6,state_mag);

for i=1:500
    for j=1:5
        d = X(i,j); state_mag = Y(i,j); z = (state_mag^2/d)^(1/2);
        state = z*ones(1,d); shift_state = ShiftState(state,d); 
        CantorCoords(i,j) = CantorPair(shift_state);
        SzudzikCoords(i,j) = SzudzikPair(shift_state);
    end
end

figure(2); clf; h=gca; set(h,'zscale','log'); hold on
surf(X,Y,CantorCoords, 'FaceColor', [1 0 0], 'DisplayName', 'Cantor Pairing');
surf(X,Y,SzudzikCoords, 'FaceColor', [0 0 1], 'DisplayName', 'Szudzik Pairing');
title('Log. Mag. of Key vs. Dim. of State vs. Mag. of State', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "northwest";
lgd.FontSize = 10;
xticks(dimension)
xlabel('Dimension of State', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('Magnitude of State', 'Interpreter', 'Latex', 'FontSize', 10)
zlabel('Log of Magnitude of Key', 'Interpreter', 'Latex', 'FontSize', 10)
view(17,22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shift_state = ShiftState(state,d)
    shift_state = zeros(1,d);
    for i=1:d
        if(state(i)<0)
            shift_state(i)=-2*state(i)-1;
        else
            shift_state(i)=2*state(i);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = CantorPair(state)
    if(length(state)>2)
        last = state(end); state(end) = [];
        x = CantorPair(state); y = last;
        key = (1/2)*(x+y)*(x+y+1)+y;
    else
        x=state(1); y=state(2);
        key = (1/2)*(x+y)*(x+y+1)+y;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = SzudzikPair(state)
    if(length(state)>2)
        last = state(end); state(end) = [];
        x = SzudzikPair(state); y = last;
        if(x<y)
            key=y^2+x;
        else
            key=x^2+x+y;
        end
    else
        x=state(1); y=state(2);
    
        if(x<y)
            key=y^2+x;
        else
            key=x^2+x+y;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%