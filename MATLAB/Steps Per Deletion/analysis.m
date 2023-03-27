file = "time_vs_steps2.txt";
blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [2 inf]);
step = blob(1,:); time = blob(2,:); 

figure(1); clf; grid on; hold on; 

title("BST Steps/Deletion Avg. Time", 'Interpreter', 'Latex');
xlabel("Steps between Deletion", 'Interpreter', 'Latex');
ylabel("Avg. Time per 2000 simualtion steps (s)", 'Interpreter', 'Latex');
scatter(step, time, 'filled'); 
yline(45.9052, '--', "Nested List Time"); 

