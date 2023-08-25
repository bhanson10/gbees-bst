x = 0:1; y = 0:1; z = 0:1;

[X, Y, Z] = meshgrid(x,y,z);
num = length(x)*length(y)*length(z);
X = reshape(X,[1,num]); Y = reshape(Y,[1,num]); Z = reshape(Z,[1,num]);

neighbors_a = {};
for i=1:num
    neighbors_a{end+1} = [X(i) Y(i) Z(i)]; 
end   

disp(neighbors_a);
disp(length(neighbors_a));

neighbors_b = {};

for i=y(1):y(end)
    for j=x(1):x(end)
        for k=z(1):z(end)
            neighbors_b{end+1} = [j i k];
        end
    end
end

disp(neighbors_b);
disp(length(neighbors_b));

disp(isequal(neighbors_a, neighbors_b));