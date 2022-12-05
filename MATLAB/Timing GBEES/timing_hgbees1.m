%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
T=1; G.thresh=0.00002; G.start=[-11.5; -10; 9.5];
G.dt=.0005; dt=.005; G.dx=0.4; G.d=3; G.sigma=4; G.b=1; G.r=48; G.L=30;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%% 
G.Y = eye(G.d,'int16'); [hD] = Initialize_D(G); 

y=G.start; ys=y; t=0; [hD]=Modify_pointset(hD,G); 

%All Plots
size = [];

%Plot 1
mod_time = []; rhs_time = [];  

%Plot 2
neighbors_time = []; vuw_time = []; remove_time = []; norm_time = [];

%Plot 3
state_conv_time = []; key_conv_time = []; key_check_time = [];

%Plot 4
initial_flux_time = []; total_flux_time = []; K_time = [];

for timestep=1:T/G.dt, disp("Timestep: " + string(timestep)); t=t+G.dt; if mod(timestep,1)==0, [hD]=Modify_pointset(hD,G); mod_t = hD.mod_t; end
    [K,hD]=RHS_P(hD,G); hD.keys = keys(hD.P); hD.P(hD.keys) = hD.P(hD.keys) + G.dt.*K; rhs_t = hD.rhs_t;
    
    %Plot 1
    mod_time = [mod_time mod_t];
    rhs_time = [rhs_time rhs_t];
    size = [size length(hD.keys)];

    %Plot 2
    neighbors_time = [neighbors_time hD.neighbors_t];
    vuw_time = [vuw_time hD.vuw_t];
    remove_time = [remove_time hD.remove_t];
    norm_time = [norm_time hD.fix_prob_t];

    %Plot 3
    state_conv_time = [state_conv_time hD.state_conv_t];
    key_conv_time = [key_conv_time hD.key_conv_t];
    key_check_time = [key_check_time hD.key_check_t];

    %Plot 4
    initial_flux_time = [initial_flux_time hD.initial_f_t];
    total_flux_time = [total_flux_time hD.total_f_t];
    K_time = [K_time hD.K_t];

end

figure(1); clf; hold on
scatter(size, mod_time, 'k', 'filled', 'DisplayName', 'Modify Pointset');
scatter(size, rhs_time, 'b', 'filled', 'DisplayName', 'RHS');
title('Contribution to Total Timestep, HGBEES1', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "best";
lgd.FontSize = 10;
xlabel('Size of Dictionary', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('time of substep (s)', 'Interpreter', 'Latex', 'FontSize', 10)

figure(2); clf; hold on
scatter(size, neighbors_time, 'k', 'filled', 'DisplayName', 'Check/Create Neighbors');
scatter(size, vuw_time, 'b', 'filled', 'DisplayName', 'Initialize VUW');
scatter(size, remove_time, 'r', 'filled', 'DisplayName', 'Remove Small Entries');
scatter(size, norm_time, 'g', 'filled', 'DisplayName', 'Normalize P');
title('Contribution to Modify Pointset, HGBEES1', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "best";
lgd.FontSize = 10;
xlabel('Size of Dictionary', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('time of substep (s)', 'Interpreter', 'Latex', 'FontSize', 10)

figure(3); clf; hold on
scatter(size, state_conv_time, 'k', 'filled', 'DisplayName', 'State Conversion');
scatter(size, key_conv_time, 'b', 'filled', 'DisplayName', 'Key Conversion');
scatter(size, key_check_time, 'r', 'filled', 'DisplayName', 'Check Key');
title('Contribution to Check/Create Neighbors, HGBEES1', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "best";
lgd.FontSize = 10;
xlabel('Size of Dictionary', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('time of substep (s)', 'Interpreter', 'Latex', 'FontSize', 10)

figure(4); clf; hold on
scatter(size, initial_flux_time, 'k', 'filled', 'DisplayName', 'Initial Flux');
scatter(size, total_flux_time, 'b', 'filled', 'DisplayName', 'Total flux');
scatter(size, K_time, 'r', 'filled', 'DisplayName', 'K');
title('Contribution to RHS, HGBEES1', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "best";
lgd.FontSize = 10;
xlabel('Size of Dictionary', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('time of substep (s)', 'Interpreter', 'Latex', 'FontSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = Initialize_D(G)
    D.P = dictionary(); D.P("0") = 0;
    D.v = dictionary(); D.v("0") = {zeros(1,G.d)};
    D.u = D.v; D.w = D.v; D.f = D.v;
   
    for i=round((G.start(1)-2)/G.dx):round((G.start(1)+2)/G.dx)
        for j=round((G.start(2)-2)/G.dx):round((G.start(2)+2)/G.dx)
            for k=round((G.start(3)-2)/G.dx):round((G.start(3)+2)/G.dx)
                state = [i j k]; key = state_conversion(state);
                x=(i*G.dx - G.start(1))^2+(j*G.dx - G.start(2))^2+(k*G.dx - G.start(3))^2;
                D.P(key) = exp(-4.*x/2.); 
            end
        end
    end
    D.keys = keys(D.P);D.values = values(D.P); 
    D.m=numEntries(D.P); D.n=D.m; D = Initialize_vuw(D,G,2); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=Initialize_vuw(D,G,b)
    for l=b:D.n
        current_key = D.keys(l);
        x=G.dx.*(key_conversion(current_key)); xh=G.dx/2;
        v1=G.sigma*(x(2)-(x(1)+xh)); v2=-(x(2)+xh)-x(1)*x(3); v3=-G.b*(x(3)+xh)+x(1)*x(2)-G.b*G.r; 
        D.v(current_key) = {[v1 v2 v3]};
        D.u(current_key)={[min(v1, 0) min(v2, 0) min(v3, 0)]};
        D.w(current_key)={[max(v1, 0) max(v2, 0) max(v3, 0)]};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(y,G)                          
    f=[G.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -G.b*y(3)+y(1)*y(2)-G.b*G.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = state_conversion(state)
    i = dec2bin(state(1),8);
    j = dec2bin(state(2),8);
    k = dec2bin(state(3),8);
    key = strcat(i,j,k);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = Modify_pointset(D,G) 
    D.keys = keys(D.P); D.m = numEntries(D.P);
    state_conv_t = 0; key_conv_t = 0; key_check_t = 0;
    for l=2:D.m    % Check/Create Neighbors of Big Cells
        if(D.P(D.keys(l))>=G.thresh)
            tic; og_state = key_conversion(D.keys(l)); key_conv_t = key_conv_t + toc;
    
            %   Faces (6)
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2), og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2), og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1), og_state(2) + 1, og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1), og_state(2) - 1, og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end   
            tic; new_key = state_conversion([og_state(1), og_state(2), og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end        
            tic; new_key = state_conversion([og_state(1), og_state(2), og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end        
     
            %   Edges (12)
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2) + 1, og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2) - 1, og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end       
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2) + 1, og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2) - 1, og_state(3)]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2), og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2), og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2), og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2), og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1), og_state(2) + 1, og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1), og_state(2) + 1, og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1), og_state(2) - 1, og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1), og_state(2) - 1, og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            
            %   Corners (8)
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2) + 1, og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2) + 1, og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2) - 1, og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2) + 1, og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) + 1, og_state(2) - 1, og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2) - 1, og_state(3) + 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2) + 1, og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
            tic; new_key = state_conversion([og_state(1) - 1, og_state(2) - 1, og_state(3) - 1]); state_conv_t = state_conv_t + toc; tic; if(~isKey(D.P, new_key)), D.P(new_key) = 0; key_check_t = key_check_t + toc; end
        end
    end
    D.state_conv_t = state_conv_t; D.key_conv_t = key_conv_t; D.key_check_t = key_check_t;
    D.neighbors_t =  state_conv_t + key_conv_t + key_check_t;

    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P);
    tic; D = Initialize_vuw(D,G,D.m+1); D.vuw_t = toc;
    tic;
    for l=2:D.m                     % Remove Small Elements
        if(D.values(l) < G.thresh)&&(no_neighbors(D,G,l))
            D.P(D.keys(l)) = []; D.v(D.keys(l)) = []; D.u(D.keys(l)) = []; D.w(D.keys(l)) = [];
        end
    end
    D.remove_t = toc;
    tic;
    D.n = numEntries(D.P); D.keys = keys(D.P);    
    D.P(D.keys) = max(D.P(D.keys), 0); D.values = values(D.P); 
    prob_sum = sum(values(D.P)); D.P(D.keys) = D.values./prob_sum; D.values = values(D.P);
    D.fix_prob_t = toc;
    D.mod_t = D.neighbors_t + D.vuw_t + D.remove_t + D.fix_prob_t;
end                                                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,D]=RHS_P(D,G)
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P); K(1:D.n,1)=0; 
    D.f(D.keys) = {zeros(1,G.d)};

    tic;
    for l=2:D.n
        l_key = D.keys(l); l_state = key_conversion(l_key);
        f_l = D.f(l_key); f_l = f_l{1}; u_l = D.u(l_key); u_l = u_l{1}; w_l = D.w(l_key); w_l = w_l{1};  
        for d=1:G.d
            k_state = l_state; k_state(d) = k_state(d)+1; k_key = state_conversion(k_state);
            if(~isKey(D.P,k_key)),k_key = "0";end
            f_l(d) = w_l(d)*D.P(l_key) + u_l(d)*D.P(k_key);
        end
        D.f(l_key) = {f_l};
    end
    D.initial_f_t = toc;

    tic;
    for d=1:G.d
        for l=2:D.n
            l_key = D.keys(l); l_state = key_conversion(l_key);
            f_l = D.f(l_key); f_l = f_l{1}; w_l = D.w(l_key); w_l = w_l{1}; 
            i_state = l_state; i_state(d) = i_state(d)-1; i_key = state_conversion(i_state);
            if(~isKey(D.P,i_key)),i_key = "0";end
            f_i = D.f(i_key); f_i = f_i{1}; w_i = D.w(i_key); w_i = w_i{1}; u_i = D.u(i_key); u_i = u_i{1};
            v_i = D.v(i_key); v_i = v_i{1};
            if (D.P(l_key)>G.thresh)||(D.P(i_key)>G.thresh)
                F=G.dt*(D.P(l_key)-D.P(i_key))/(2*G.dx);
                for e=1:G.d
                    if e~=d
                        j_state = l_state; j_state(e) = j_state(e)-1; j_key = state_conversion(j_state);
                        if(~isKey(D.P,j_key)),j_key = "0";end
                        f_j = D.f(j_key); f_j = f_j{1}; w_j = D.w(j_key); w_j = w_j{1}; u_j = D.u(j_key); u_j = u_j{1};
                        p_state = i_state; p_state(e) = p_state(e)-1; p_key = state_conversion(p_state);
                        if(~isKey(D.P,p_key)),p_key = "0";end
                        f_p = D.f(p_key); f_p = f_p{1}; w_p = D.w(p_key); w_p = w_p{1}; u_p = D.u(p_key); u_p = u_p{1};

                        f_l(e) = f_l(e)-w_l(e)*w_i(d)*F;
                        f_j(e) = f_j(e)-u_j(e)*w_i(d)*F; D.f(j_key) = {f_j};
                        f_i(e) = f_i(e)-w_i(e)*u_i(d)*F; 
                        f_p(e) = f_p(e)-u_p(e)*u_i(d)*F; D.f(p_key) = {f_p};
                    end
                end
                D.f(l_key) = {f_l};
                D.f(i_key) = {f_i};
                
                i_i_state = i_state; i_i_state(d) = i_i_state(d)-1; i_i_key = state_conversion(i_i_state);
                if(~isKey(D.P,i_i_key)),i_i_key = "0";end
                k_state = l_state; k_state(d) = k_state(d)+1; k_key = state_conversion(k_state);
                if(~isKey(D.P,k_key)),k_key = "0";end
                if(v_i(d)>0)
                    th=(D.P(i_key)-D.P(i_i_key))/(D.P(l_key)-D.P(i_key));
                else
                    th=(D.P(k_key)-D.P(l_key))/(D.P(l_key)-D.P(i_key));
                end
                
                f_i = D.f(i_key); f_i = f_i{1};
                t=abs(v_i(d)); f_i(d)=f_i(d)+t*(G.dx/G.dt-t)*F*MC(th); D.f(i_key) = {f_i};
            end
        end
    end
    D.total_f_t = toc;
    
    tic;
    for l=2:D.n
        l_key = D.keys(l); l_state = key_conversion(l_key);
        f_l = D.f(l_key); f_l = f_l{1};
        for d=1:G.d
            i_state = l_state; i_state(d) = i_state(d)-1; i_key = state_conversion(i_state);
            if(~isKey(D.P,i_key)),i_key = "0";end
            f_i = D.f(i_key); f_i = f_i{1};
            K(l,1)=K(l,1)-(f_l(d)-f_i(d))/G.dx;
        end
    end
    D.K_t = toc;
    D.rhs_t = D.initial_f_t + D.total_f_t + D.K_t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi]=MC(th), phi=max(0,min([(1+th)/2 2 2*th]));   end   % Flux limiters
function [phi]=VL(th), phi=min((th+abs(th))/(1+abs(th)),0); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool=no_neighbors(D,G,l)
    bool = 1;
    current_key = D.keys(l); og_state = key_conversion(current_key);
    % Faces (6)
    new_key = state_conversion([og_state(1)+1 og_state(2) og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2) og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)+1 og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)-1 og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = state_conversion([og_state(1) og_state(2) og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = state_conversion([og_state(1) og_state(2) og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    
    % Edges (12)
    new_key = state_conversion([og_state(1)+1 og_state(2)+1 og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)-1 og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = state_conversion([og_state(1)-1 og_state(2)+1 og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end   
    new_key = state_conversion([og_state(1)-1 og_state(2)-1 og_state(3)]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2) og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2) og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2) og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2) og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)+1 og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)+1 og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)-1 og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)-1 og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    %   Corners (8)
    new_key = state_conversion([og_state(1)+1 og_state(2)+1 og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)+1 og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)-1 og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)+1 og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)-1 og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)-1 og_state(3)+1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)+1 og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)-1 og_state(3)-1]); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
end