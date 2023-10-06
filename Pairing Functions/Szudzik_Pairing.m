state = [1 1 1 1 1 1], d=length(state);
shift_state = ShiftState(state,d),
key = SzudzikPair(shift_state),
shift_state = SzudzikUnpair(key,d),
state = UnshiftState(shift_state),

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
function key = SzudzikPair(state)
    if(length(state)>2)
        last = state(end); state(end) = [];
        x = SzudzikPair(state); y = last;
        if(x<y)
            key=y^2+x;
        else
            key=x^2+x+y;
        end
    else
        x=state(1); y=state(2);
    
        if(x<y)
            key=y^2+x;
        else
            key=x^2+x+y;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = SzudzikUnpair(key,d)
    state = zeros(1,d);
    for i=2:d
        z=key; m=floor(z^(1/2));
        if((z-m^2)<m)
            state(d-i+1) = z-m^2; state(d-i+2) = m;
        else
            state(d-i+1) = m; state(d-i+2) = z-m^2-m;
        end
        key = state(d-i+1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = UnshiftState(shift_state)
    d = length(shift_state); state = zeros(1,d);
    for i=1:d
        if(mod(shift_state(i),2)==0)
            state(i)=shift_state(i)/2;
        else
            state(i)=(shift_state(i)+1)/-2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%