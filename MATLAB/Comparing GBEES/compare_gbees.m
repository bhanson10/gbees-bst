function compare_gbees(D, hD)
    d_states = {}; d_prob = []; d_v = {}; d_u = {}; d_w = {}; d_f = {};
    f_initialized = any(D.f(:));
    for l=2:D.m
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
    
    v1 = hD.v{1}; v2 = hD.v{2}; v3 = hD.v{3};
    u1 = hD.u{1}; u2 = hD.u{2}; u3 = hD.u{3}; 
    w1 = hD.w{1}; w2 = hD.w{2}; w3 = hD.w{3}; 
    f1 = hD.f{1}; f2 = hD.f{2}; f3 = hD.f{3}; 
    hf_initialized = ~(numEntries(f1) == 0);

    for l=1:hD.n
        hd_states{end+1} = key_conversion(hD.keys(l));
        hd_prob = [hd_prob hD.P(hD.keys(l))]; 
        hd_v{end+1} = [v1(hD.keys(l)) v2(hD.keys(l)) v3(hD.keys(l))];
        hd_u{end+1} = [u1(hD.keys(l)) u2(hD.keys(l)) u3(hD.keys(l))];
        hd_w{end+1} = [w1(hD.keys(l)) w2(hD.keys(l)) w3(hD.keys(l))];
        if (hf_initialized)
            hd_f{end+1} = [f1(hD.keys(l)) f2(hD.keys(l)) f3(hD.keys(l))];
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