hash_table = dictionary();
for i=1:5
    state = round(-5 + 10*rand(1, 6), 1);
    key = strjoin(string(state), '\');
    value = rand(1);
    hash_table(key) = value;
end

disp(hash_table);