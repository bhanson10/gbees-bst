hash_table = dictionary();
for i=1:5
    key = str2num(strjoin(string(randi([0 1],1,5)), ''));
    value = rand(1);
    hash_table(key) = value;
end

disp(hash_table);