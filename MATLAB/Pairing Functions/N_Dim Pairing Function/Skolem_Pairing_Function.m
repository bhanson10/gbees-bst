%%%%%%%%%%%%%%%%%%%% Examining every state on a Shell %%%%%%%%%%%%%%%%%%%%%
%{
s = 1; n = 3; % Inputs
num_states = f(s, n); 
states = g(s,n);
shift_states = Shift(states, n);
keys = skolem(shift_states);
[keys,sortIdx] = sort(keys,'ascend');
states = states(sortIdx);
print_results(states, keys, n);
%}
%%%%%%%%%%%%%%%%%%%% Examining every state on a Shell %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Getting from state to key and back %%%%%%%%%%%%%%%%%%%%
states = {[10 10 10 10 10 10]}; n = length(states{1});
shift_states = Shift(states, n);
keys = skolem(shift_states);
print_results(states, keys, n);
%%%%%%%%%%%%%%%%%%% Getting from state to key and back %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scaled_vals = scale(keys, num_states)
    scaled_vals = zeros(1,length(keys));
    r_max = keys(end); r_min = keys(1);
    t_max = num_states; t_min = 1; 

    for i=1:length(keys)
        scaled_vals(1,i) = (keys(i)-r_min)*(t_max-t_min)/(r_max-r_min)+t_min;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function num = f(s,n)
    if(s == 0)
        num = 1;
    else
        sum_prev_s = 0;
        for i=0:s-1
            sum_prev_s = sum_prev_s + f(i,n);
        end
        num = (2*s+1)^n - sum_prev_s; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shift_states = Shift(states,n)
    shift_states = {}; 
    for i=1:length(states)
        state = states{i};
        shift_state = zeros(1,n);
        for j=1:n
            if state(j) >= 0
                shift_state(1,j) = 2*state(j)+1;
                
            else
                shift_state(1,j) = (-2*state(j));
            end
        end
        shift_states{end+1} = shift_state; 
    end
end           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function skolem_vals = skolem(shift_states)
    skolem_vals = zeros(1,length(shift_states));

    for i=1:length(shift_states)
        state = shift_states{i};
        
        skolem_val = 0;
        for j=1:length(state)
            total = sum(state(1:j)); 
            skolem_val = skolem_val + nchoosek(total+j-1, j);
        end

        skolem_vals(1,i) = skolem_val; 
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
function print_results(states, keys, n)

    switch n
        case 1
            for i=1:length(states)
                state = states{i};
                fprintf("State: (%d)    Key: %d\n", state(1), keys(i));
            end
        case 2
            for i=1:length(states)
                state = states{i};
                fprintf("State: (%d, %d)    Key: %d\n", state(1), state(2), keys(i));
            end
        case 3
            for i=1:length(states)
                state = states{i};
                fprintf("State: (%d, %d, %d)    Key: %d\n", state(1), state(2), state(3), keys(i));
            end
        case 4
            for i=1:length(states)
                state = states{i};
                fprintf("State: (%d, %d, %d, %d)    Key: %d\n", state(1), state(2), state(3), state(4), keys(i));
            end
        case 5
            for i=1:length(states)
                state = states{i};
                fprintf("State: (%d, %d, %d, %d, %d)    Key: %d\n", state(1), state(2), state(3), state(4), state(5), keys(i));
            end
        case 6
            for i=1:length(states)
                state = states{i};
                fprintf("State: (%d, %d, %d, %d, %d, %d)    Key: %d\n", state(1), state(2), state(3), state(4), state(5), state(6), keys(i));
            end
        otherwise
            disp("Inapplicable dimension.")
    end
end