file = "time_complexity.txt";
blob_ID = fopen(file, 'r'); blob = fscanf(blob_ID, '%f', [4 inf]);
n = blob(1,:); t = blob(2,:)'; 

figure(1); hold on
title("Time Complexity of HGBEES")
xlabel("Hash Table Size")
ylabel("Time (s)")
scatter(n,t)
