clearvars
state = [0 0]; d = length(state); v = [1 1]; count = 1; 

%{
for a=1:d
    x_state = state;
    if(v(a) > 0)
         x_state(a) = x_state(a) + 1; neighbors(1,count) = state2str(x_state); count = count + 1; 
         for b=1:d-1
             if(b~=a)
                 if(v(b)>0)
                     y_state = x_state; y_state(b) = y_state(b) + 1; 
                     neighbors(1,count) = state2str(y_state); count = count + 1; 
                     for c=1:d-2
                         if((c~=a)&&(c~=b))
                             if(v(c)>0)
                                 z_state = y_state; z_state(c) = z_state(c) + 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             else
                                 z_state = y_state; z_state(c) = z_state(c) - 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             end
                         end
                     end
                 else
                     y_state = x_state; y_state(b) = y_state(b) - 1; 
                     neighbors(1,count) = state2str(y_state); count = count + 1; 
                     for c=1:d-2
                         if((c~=a)&&(c~=b))
                             if(v(c)>0)
                                 z_state = y_state; z_state(c) = z_state(c) - 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             else
                                 z_state = y_state; z_state(c) = z_state(c) - 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             end
                         end
                     end
                 end
             end
         end
    else
        x_state(a) = x_state(a) - 1; neighbors(1,count) = state2str(x_state); count = count + 1; 
        for b=1:d-1
             if(b~=a)
                 if(v(b)>0)
                     y_state = x_state; y_state(b) = y_state(b) + 1; 
                     neighbors(1,count) = state2str(y_state); count = count + 1; 
                     for c=1:d-2
                         if((c~=a)&&(c~=b))
                             if(v(c)>0)
                                 z_state = y_state; z_state(c) = z_state(c) + 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             else
                                 z_state = y_state; z_state(c) = z_state(c) - 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             end
                         end
                     end
                 else
                     y_state = x_state; y_state(b) = y_state(b) - 1; 
                     neighbors(1,count) = state2str(y_state); count = count + 1; 
                     for c=1:d-1
                         if((c~=a)&&(c~=b))
                             if(v(c)>0)
                                 z_state = y_state; z_state(c) = z_state(c) + 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             else
                                 z_state = y_state; z_state(c) = z_state(c) - 1;
                                 neighbors(1,count) = state2str(z_state); count = count + 1; 
                             end
                         end
                     end
                 end
             end
         end
    end
end

disp(neighbors);
%}

iter = 2; count = 1; anti_iter = 1; 
recursive_for(state, d, iter, anti_iter, v);
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
function recursive_for(state,d,iter,anti_iter,v)
    if (iter == 0)
        return 
    else
        for a=1:iter
            new_state = state; 
            if(v(anti_iter)>0)
                new_state(anti_iter) = new_state(anti_iter) + 1; disp(state2str(new_state));  
                recursive_for(new_state, d, iter-1, anti_iter+1, v)
            else
                new_state(anti_iter) = new_state(anti_iter) - 1; disp(state2str(new_state)); 
                recursive_for(new_state, d, iter-1, anti_iter+1, v)
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%