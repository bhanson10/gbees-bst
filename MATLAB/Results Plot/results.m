time = [0 400 800 1200 1600 2000];
nl = [0 6.822 8.248 11.638 20.103 43.983];
ht = [0 24.576 52.547 89.345 180.273 409.829];
bst = [0 11.774 25.789 43.404 103.885 302.746];

figure(1); hold on; grid on; 
xlabel("Timestep");
ylabel("Total time (s)");
plot(time, nl, 'r', 'DisplayName', 'Nested Lists');
plot(time, ht, 'b', 'DisplayName', 'Hash Table');
plot(time, bst, 'g', 'DisplayName', 'BST');
lgd = legend;
lgd.Location = "northwest";
lgd.FontSize = 10;