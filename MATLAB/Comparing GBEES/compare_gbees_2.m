function compare_gbees_2(D, hD)
    d_states = {}; d_prob = []; d_v = {}; d_u = {}; d_w = {}; d_f = {};
    f_initialized = any(D.f(:));
    for l=2:D.n
        d_states{end+1} = [D.j(l,1) D.j(l,2) D.j(l,3)];
        d_prob = [d_prob D.P(l)];
        d_v{end+1} = [D.v(l,1) D.v(l,2) D.v(l,3)];
        d_u{end+1} = [D.u(l,1) D.u(l,2) D.u(l,3)];
        d_w{end+1} = [D.w(l,1) D.w(l,2) D.w(l,3)];
        if (f_initialized)
            d_f{end+1} = [D.f(l,1) D.f(l,2) D.f(l,3)];
        end
    end

    hd_states = {}; hd_prob = []; hd_v = {}; hd_u = {}; hd_w = {}; hd_f = {};
    hf_initialized = ~(numEntries(hD.f) == 1);

    for l=2:hD.n
        state = hD.keys(l); state = state{1};
        hd_states{end+1} = state;
        hd_prob = [hd_prob hD.P(hD.keys(l))]; 
        v = hD.v(hD.keys(l)); v = v{1};
        hd_v{end+1} = v;
        u = hD.u(hD.keys(l)); u = u{1};
        hd_u{end+1} = u;
        w = hD.w(hD.keys(l)); w = w{1};
        hd_w{end+1} = w;
        if (hf_initialized)
            f = hD.f(hD.keys(l)); f = f{1};
            hd_f{end+1} = f;
        end
    end
    
    match = isequal(d_states, hd_states);
    if(match)
        disp("States match.")
    else
        disp("States do not match.")
        disp("# of D states: " + string(length(d_states)));
        disp("# of hD states: " + string(length(hd_states)));
    end

    match = isequal(d_prob, hd_prob);
    if(match)
        disp("Probabilities match.")
    else
        disp("Probabilities do not match.")
    end

    match = isequal(d_v, hd_v);
    if(match)
        disp("v matches.")
    else
        disp("v do not match.")
    end

    match = isequal(d_u, hd_u);
    if(match)
        disp("u matches.")
    else
        disp("u do not match.")
    end

    match = isequal(d_w, hd_w);
    if(match)
        disp("w matches.")
    else
        disp("w do not match.")
    end
    
    if(f_initialized)&&(hf_initialized)
        match = isequal(d_f, hd_f);
        if(match)
            disp("f matches.")
        else
            disp("f do not match.")
        end
    end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = key_conversion(input_key)
    input_key = convertStringsToChars(input_key);
    length = 8;
    bin_i = input_key(1:8); bin_j = input_key(9:16); bin_k = input_key(17:24);
    i = twos2dec(bin_i, length); j = twos2dec(bin_j, length); k = twos2dec(bin_k, length);
    state = [i j k];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decimal = twos2dec(x,bits)
    if (x(1) == '0')
        decimal = bin2dec(x);
    else
        for i=1:bits
            if(x(i) == '0') 
                x(i) = '1';
            elseif(x(i) == '1')
                x(i) = '0';
            end
        end
        decimal = -bin2dec(dec2bin(bin2dec(x) + bin2dec('1')));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%