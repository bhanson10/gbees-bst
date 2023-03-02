file1 = "bst.txt"; file2 = "ht.txt"; G.d = 3;

blob_ID = fopen(file1, 'r'); blob = fscanf(blob_ID, '%f', [4 inf]);
bst_probs = blob(1,:); bst_set = blob(2:4,:)'; bst_n = length(bst_probs);

blob_ID = fopen(file2, 'r'); blob = fscanf(blob_ID, '%f', [4 inf]);
ht_probs = blob(1,:); ht_set = blob(2:4,:)'; ht_n = length(ht_probs);

bst_keys = zeros(1,bst_n); bst_state = zeros(1,3);

for i=1:bst_n
    bst_state(1) = bst_set(i,1); 
    bst_state(2) = bst_set(i,2);
    bst_state(3) = bst_set(i,3);
    key = state_conversion(bst_state,G);
    bst_keys(i) = key; 
end

ht_keys = zeros(1,ht_n); ht_state = zeros(1,3);

for i=1:ht_n
    ht_state(1) = ht_set(i,1); 
    ht_state(2) = ht_set(i,2);
    ht_state(3) = ht_set(i,3);
    key = state_conversion(ht_state,G);
    ht_keys(i) = key; 
end

C = setxor(bst_keys, ht_keys); 
disp(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = state_conversion(state,G)
    shift_state = ShiftState(state,G.d);
    key = CantorPair(shift_state);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = key_conversion(key,G)
    shift_state = CantorUnpair(key,G.d);
    state = UnshiftState(shift_state);
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
function key = CantorPair(state)
    if(length(state)>2)
        last = state(end); state(end) = [];
        x = CantorPair(state); y = last;
        key = (1/2)*(x+y)*(x+y+1)+y;
    else
        x=state(1); y=state(2);
        key = (1/2)*(x+y)*(x+y+1)+y;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = CantorUnpair(key,d)
    state = zeros(1,d);
    for i=2:d
        z=key; w=floor(((8*z+1)^(1/2)-1)/2); t=(w^2+w)/2; 
        y=z-t; x=w-y; state(d-i+1)=x; state(d-i+2) =y;
        key=x;
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