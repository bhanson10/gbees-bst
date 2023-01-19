d = 10000; symbols = ['a':'z' 'A':'Z']; st_length = zeros(1,d); 
create_times = zeros(1,d); search_times = zeros(1,d);

for i=1:d % Testing Dictonaries of Length:1-1000
    dic = dictionary(); string_length = i; st_length(i) = string_length;
    create_time = 0;
    for j=1:100 % Creating 100 entries
        key = symbols(randi(numel(symbols),[1 string_length])); value = j;
        tic;
        dic(key) = value;
        total = toc;
        create_time = create_time + total;
    end
    create_times(i) = create_time;
    tic;
    for j=1:100 %Performing 100 Searches
        key = string(j);
        isKey(dic,key);
    end
    search_times(i) = toc;
end

figure(1); clf; hold on, 
title("Creating a 100 entry dictionary, searching 100 times")
plot(st_length, create_times, 'r', "DisplayName", "Creating Dictionary")
plot(st_length, search_times, 'b', "DisplayName", "Searching Dictionary")
xlabel("Key Length")
ylabel("Time (s)")
lgd = legend;
lgd.Location = "best";
lgd.FontSize = 10;
