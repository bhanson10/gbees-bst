max = 10;
counts = [];
e_times = [];
ne_times = [];
for i = 1:max
    hash_table = dictionary();
    range = i;
    count = 0;
    for pos_x = 1:range
        for pos_y = 1:range
            for pos_z = 1:range
                for vel_x = 1:range
                    for vel_y = 1:range
                        for vel_z = 1:range
                            count = count + 1;
                            state =  [pos_x pos_y pos_z 
                                vel_x vel_y vel_z];
                            key = hash(state);
                            prob = rand(1);
                            hash_table(key) = prob;
                        end
                    end
                end
            end
        end
    end
    
    counts = [counts count];

    % For Existing Keys
    e_rand_state = randi([1, range], 1,6);
    e_rand_key = hash(e_rand_state);
   
    total_time = 0;
    for j=1:10
        tic;
        isKey(hash_table,e_rand_key);
        elapsed_time = toc;
        total_time = total_time + elapsed_time;
    end
    avg_time = total_time/10;
    e_times = [e_times (avg_time*1000)];

    % For Non-Existing Keys
    ne_rand_state = randi([range + 1, range + 10], 1,6);
    ne_rand_key = hash(ne_rand_state);
   
    total_time = 0;
    for j=1:10
        tic;
        isKey(hash_table,ne_rand_key);
        elapsed_time = toc;
        total_time = total_time + elapsed_time;
    end
    avg_time = total_time/10;
    ne_times = [ne_times (avg_time*1000)];
end

figure; clf;
set(gca,'fontsize', 18); 
semilogx(counts,e_times,'o','MarkerFaceColor',[0 0.447 0.741], 'DisplayName', 'Existing States')
hold on;
semilogx(counts,ne_times,'o','MarkerFaceColor',[1 0.5 0], 'DisplayName','Non-existing States')
grid on
title("Hash Table Search for State");
xlabel("# of Hash Table Elements")
ylabel("Time to search (ms)")
legend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------Functions---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = hash(state);
    for i=1:length(state)
        if i == 1
            key = string(state(i));
        else
            key = key + '/' + string(state(i));
        end
    end
end 



                  