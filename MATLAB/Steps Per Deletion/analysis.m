file = "time_vs_per.txt";
blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [2 inf]);
step = blob(1,:); time = blob(2,:); 

figure(1); clf; grid on; hold on; 
xlim([0,100])
title("BST \%/Deletion Avg. Time", 'Interpreter', 'Latex');
xlabel("\% of Inactive Cells before Deletion", 'Interpreter', 'Latex');
ylabel("Avg. Time per 2000 simualtion steps (s)", 'Interpreter', 'Latex');
plot(step,time,'--blue');
scatter(step, time, 'blue','filled'); 
yline(33.77, '--', "Nested List Time"); 

