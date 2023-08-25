function compare_gbees(hD1, hD2)
    hd1_states = {}; hd1_prob = []; hd1_v = {}; hd1_u = {}; hd1_w = {}; hd1_f = {};
    
    v = hD1.v; u = hD1.u; w = hD1.w; f = hD1.f; 
    hf_initialized = ~(numEntries(f) == 1); 

    for l=2:hD1.n
        hd1_states{end+1} = key_conversion1(hD1.keys(l));
        hd1_prob = [hd1_prob hD1.P(hD1.keys(l))]; 
        current_v = v(hD1.keys(l)); current_v = current_v{1};
        hd1_v{end+1} = [current_v(1) current_v(2) current_v(3)];
        current_u = u(hD1.keys(l)); current_u = current_u{1};
        hd1_u{end+1} = [current_u(1) current_u(2) current_u(3)];
        current_w = w(hD1.keys(l)); current_w = current_w{1};
        hd1_w{end+1} = [current_w(1) current_w(2) current_w(3)];

        if (hf_initialized)
            current_f = f(hD1.keys(l)); current_f = current_f{1};
            hd1_f{end+1} = [current_f(1) current_f(2) current_f(3)];
        end
    end
    
    G.dt=.0005; dt=.005; G.dx=0.4; G.d=3; G.sigma=4; G.b=1; G.r=48; G.L=30;
    G.N_bits = 8; G.N_data = G.d; G.fac=uint64(2^G.N_bits); G.offset32=int32(G.fac/2); G.offset64=int64(G.offset32);
    hd2_states = {}; hd2_prob = []; hd2_v = {}; hd2_u = {}; hd2_w = {}; hd2_f = {};
    
    v = hD2.v; u = hD2.u; w = hD2.w; f = hD2.f;
    hf_initialized = ~(numEntries(f) == 1);

    for l=2:hD2.n
        hd2_states{end+1} = key_conversion2(hD2.keys(l),G);
        hd2_prob = [hd2_prob hD2.P(hD2.keys(l))]; 
        current_v = v(hD2.keys(l)); current_v = current_v{1};
        hd2_v{end+1} = [current_v(1) current_v(2) current_v(3)];
        current_u = u(hD2.keys(l)); current_u = current_u{1};
        hd2_u{end+1} = [current_u(1) current_u(2) current_u(3)];
        current_w = w(hD2.keys(l)); current_w = current_w{1};
        hd2_w{end+1} = [current_w(1) current_w(2) current_w(3)];

        if (hf_initialized)
            current_f = f(hD2.keys(l)); current_f = current_f{1};
            hd2_f{end+1} = [current_f(1) current_f(2) current_f(3)];
        end
    end

    
    match = isequal(hd1_states, hd2_states);
    if(match)
        disp("States match.")
    else
        disp("States do not match.")
        disp("# of hD1 states: " + string(length(hd1_states)));
        disp("# of hD2 states: " + string(length(hd2_states)));
    end

    match = isequal(hd1_prob, hd2_prob);
    if(match)
        disp("Probabilities match.")
    else
        disp("Probabilities do not match.")
    end

    match = isequal(hd1_v, hd2_v);
    if(match)
        disp("v matches.")
    else
        disp("v do not match.")
    end

    match = isequal(hd1_u, hd2_u);
    if(match)
        disp("u matches.")
    else
        disp("u do not match.")
    end

    match = isequal(hd1_w, hd2_w);
    if(match)
        disp("w matches.")
    else
        disp("w do not match.")
    end
    
    if(f_initialized)&&(hf_initialized)
        match = isequal(hd1_f, hd2_f);
        if(match)
            disp("f matches.")
        else
            disp("f do not match.")
        end
    end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = key_conversion1(input_key)
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
function state = key_conversion2(key,G)
    for i=G.N_data:-1:1
        state(i)=idivide(key,G.fac^(i-1),'floor');
        key=key-state(i)*G.fac^(i-1);
    end
    state=int64(state)-G.offset64;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%