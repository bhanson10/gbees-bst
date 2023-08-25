%%%%%%%%%%%%%%%%%%%% Examining every state on a Shell %%%%%%%%%%%%%%%%%%%%%
%{
s = 2; n = 3; % Inputs
states = g(s,n); num = length(states);
keys = zeros(1,num);
for i=1:num
    state = states{i};
    keys(1,i) = state_conversion(state,n,0);
end
[keys,sortIdx] = sort(keys,'ascend');
states = states(sortIdx);
print_results(states, keys, n);
%}
%%%%%%%%%%%%%%%%%%%% Examining every state on a Shell %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Examining Individual States %%%%%%%%%%%%%%%%%%%%%%%%
format short; 
state = [0 109 0 0], % Inputs
d = length(state);
format long; 
key = state_conversion(state,d,1),
format short; 
state = key_conversion(key,d,1),
%%%%%%%%%%%%%%%%%%%%%% Examining Individual States %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = state_conversion(state,d,shift_status)
    if shift_status == 1
        shift_state = ShiftState(state,d),
    else
        shift_state = state;
    end
    m = max(shift_state),
    key = RosePair(shift_state,m);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = key_conversion(key,d,shift_status)
    m = floor(key^(1/d));
    shift_state = RoseUnpair(key,m,d),
    if shift_status == 1
        state = UnshiftState(shift_state,d);
    else
        state = shift_state;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shift_state = ShiftState(state,d)
    shift_state = zeros(1,d);
    for i=1:d
        if(state(i)<0)
            shift_state(i)=-2*state(i)-1;
        else
            shift_state(i)=2*state(i);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = RosePair(state,m)
    d = length(state);
    if(length(state)==1)
        key = state(1);
    else
        new_m = max(state(1:d-1));
        key = RosePair(state(1:d-1),new_m) + m^d + (m-state(d))*((m+1)^(d-1)-m^(d-1)); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = RoseUnpair(key,m,d)
    if d == 1
        state = [key];
    else
        x_d = m - floor(max(0,key-m^d-m^(d-1))/((m+1)^(d-1)-m^(d-1))); 
        new_key = key-(m^d)-(m-x_d)*((m+1)^(d-1)-m^(d-1));
        new_m = floor(new_key^(1/(d-1))); 
        state = [RoseUnpair(new_key, new_m, d-1) x_d];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = UnshiftState(shift_state,d)
    state = zeros(1,d);
    for i=1:d
        if(mod(shift_state(i),2)==0)
            state(i)=shift_state(i)/2;
        else
            state(i)=(shift_state(i)+1)/-2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function states = g(s,n)
    pos = [0:s]; states = {};

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