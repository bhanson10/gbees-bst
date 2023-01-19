%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 2; n = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_states = f(s, n); 
states = g(s,n);
keys = f(s-1,n)+[1:num_states];

for i=1:num_states
    fprintf("State: ");
    fprintf('%g', states{i});
    fprintf(" ");
    fprintf("Key: ");
    fprintf('%d', keys(i));
    fprintf('\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = f(s,n)
    if(s == 0)
        n = 1;
    else
        sum_prev_s = 0;
        for i=0:s-1
            sum_prev_s = sum_prev_s + f(i,n);
        end
        n = (2*s+1)^n - sum_prev_s; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function states = g(s,n)
    pos = [-s:s]; num = 2*s+1; states = {};

    switch n
        case 1
            for i=pos
                state = i; 
                if (max(abs(state)) == s)
                    states{end+1} = state
                end
            end
        case 2 
            for i=pos
                for j=pos 
                    state = [i j];
                    if (max(abs(state)) == s)
                        states{end+1} = state;
                    end
                end
            end
        case 3
            for i=pos
                for j=pos 
                    for k=pos
                        state = [i j k];
                        if (max(abs(state)) == s)
                            states{end+1} = state;
                        end
                    end
                end
            end
        case 4
            for i=pos
                for j=pos 
                    for k=pos
                        for l=pos
                            state = [i j k l];
                            if (max(abs(state)) == s)
                                states{end+1} = state;
                            end
                        end
                    end
                end
            end
        case 5
            for i=pos
                for j=pos 
                    for k=pos
                        for l=pos
                            for m=pos
                                state = [i j k l m];
                                if (max(abs(state)) == s)
                                    states{end+1} = state;
                                end
                            end
                        end
                    end
                end
            end
        case 6
            for i=pos
                for j=pos 
                    for k=pos
                        for l=pos
                            for m=pos
                                for o=pos
                                    state = [i j k l m o];
                                    if (max(abs(state)) == s)
                                        states{end+1} = state;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        otherwise
            disp("Inapplicable dimension.")
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
