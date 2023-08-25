file = "time_complexity_bst.txt";
blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [2 inf]);
step = blob(1,:); time = blob(2,:); 

file2 = "time_complexity_nl.txt";
blob_ID2 = fopen(file2, 'r'); blob2 = fscanf(blob_ID2, '%f', [2 inf]);
step2 = blob2(1,:); time2 = blob2(2,:); 

figure(1); clf; grid on; hold on; 
title("BST Time Complexity", 'Interpreter', 'Latex');
xlabel("Number of Cells", 'Interpreter', 'Latex');
ylabel("Time of step (s)", 'Interpreter', 'Latex');
scatter(step, time, 'blue','filled'); 
scatter(step2, time2, 'red','filled'); 