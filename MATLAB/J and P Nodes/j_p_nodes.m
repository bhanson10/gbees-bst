clearvars
state = [0 0 0]; d = length(state); 

for q=1:d
    i_state = state; i_state(q) = i_state(q)-1; i_nodes(1,q) = state2str(i_state);
    for e=1:d
        if(e~=q)
            j_nodes(e,q) = state2str(i_state); 
            if(e>q)
                p_state = i_state; p_state(e) = p_state(e) - 1;
                node = state2str(p_state);
                p_nodes(q,e) = node;
            end   
        end
    end
end

% Make p matrix diagonally symmetric
for q=1:d
    for e=1:d
        if(e<q)
            p_nodes(q,e) = p_nodes(e,q);
        end
    end
end

disp("i nodes"); 
disp(i_nodes)
disp(" ");
disp("j nodes"); 
disp(j_nodes)
disp(" ");
disp("p nodes"); 
disp(p_nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function string = state2str(state)
    string = "(";
    for i=1:length(state)
        if i == length(state)
            string = string + num2str(state(i)) + ")";
        else
            string = string + num2str(state(i)) + ", ";
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%