clearvars
state = [0 0 0]; v = [1 1 1]; d = length(state); 

iter = d; non = zeros(1,d); 
get_neighbors(state, iter, v, non);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_neighbors(state, iter, v, non)
    if (iter == 0)
        return
    else
        for a=1:iter
            if (~ismember(a,non))
                non(1,iter) = a; 
                if(v(a)>0)
                    new_state = state; new_state(a) = new_state(a) + 1; 
                    disp(state2str(new_state)); 
                    get_neighbors(new_state, iter-1, v, non);
                else
                    new_state = state; new_state(a) = new_state(a) - 1;
                    disp(state2str(new_state));
                    get_neighbors(new_state, iter-1, v, non); 
                end
            end
        end
    end   
end
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