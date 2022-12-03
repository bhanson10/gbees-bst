hD = hgbees;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hD = hgbees
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
T=1; G.thresh=0.00002; G.start=[-11.5; -10; 9.5];
G.dt=.0005; dt=.005; G.dx=0.4; G.d=3; G.sigma=4; G.b=1; G.r=48; G.L=30;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%% 
G.Y = eye(G.d,'int16'); [hD] = Initialize_D(G); 

y=G.start; ys=y; t=0;disp("First Modify pointset"); [hD]=Modify_pointset(hD,G); 

for timestep=1, t=t+G.dt; if mod(timestep,1)==0, disp("Second Modify pointset"); [hD]=Modify_pointset(hD,G); end
  K=RHS_P(hD,G); tic; hD.values = values(hD.P); hD.keys = keys(hD.P); hD.P(hD.keys) = hD.values + G.dt.*K; disp("Updating Probability: " + string(toc));
  tic; k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
  ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew; disp("Time-marching y: " + string(toc));
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = Initialize_D(G)
    D.P = dictionary();

    tic;
    D.v = {dictionary(),dictionary(),dictionary()};
    D.u = {dictionary(),dictionary(),dictionary()};
    D.w = {dictionary(),dictionary(),dictionary()};
    D.f = {dictionary(),dictionary(),dictionary()};
    disp("Initialize D - Initializing Variables: " + string(toc));

    tic;
    for i=round((G.start(1)-2)/G.dx):round((G.start(1)+2)/G.dx)
        for j=round((G.start(2)-2)/G.dx):round((G.start(2)+2)/G.dx)
            for k=round((G.start(3)-2)/G.dx):round((G.start(3)+2)/G.dx)
                key = {[i j k]};
                x=(i*G.dx - G.start(1))^2+(j*G.dx - G.start(2))^2+(k*G.dx - G.start(3))^2;
                D.P(key) = exp(-4.*x/2.); 
            end
        end
    end
    disp("Setting initial dict. entries: " + string(toc));

    D.keys = keys(D.P);D.values = values(D.P); 
    D.m=numEntries(D.P); D.n=D.m; D = Update_vuw(D,G,1); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=Update_vuw(D,G,b)
    v_time = 0; u_time = 0; w_time = 0;
    for l=b:D.n
        current_key = D.keys(l);
        x=G.dx.*(current_key{1}); xh=G.dx/2;
        tic;
        v1 = D.v{1}; v2 = D.v{2}; v3 = D.v{3};
        v1(current_key) = G.sigma*(x(2)-(x(1)+xh));
        v2(current_key) = -(x(2)+xh)-x(1)*x(3); 
        v3(current_key) = -G.b*(x(3)+xh)+x(1)*x(2)-G.b*G.r;
        v_time = v_time + toc;
        D.v = {v1, v2, v3};
        for d=1:G.d
            % ----------------------------- 3D Solid Body Rotation ------------------------------
            %D.v(l,1)=2.*D.j(2)*G.dx; D.v(l,2)=-2.*D.j(1)*G.dx; D.v(l,3)=-.5;
            % ----------------------------------- 3D Lorenz -------------------------------------
            v = D.v{d}; w = D.w{d}; u = D.u{d}; 
            tic;
            w(current_key)=max(v(current_key), 0); w_time = w_time + toc;
            tic;
            u(current_key) = min(v(current_key), 0); u_time = u_time + toc;
            D.u{d} = u; D.w{d} = w;
            % ------------------------- (keep one of the above 2 sections) ----------------------
        end
    end
    disp("Update_vuw - v time: " + string(v_time));
    disp("Update_vuw - u time: " + string(u_time));
    disp("Update_vuw - w time: " + string(w_time));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(y,G)                          
    f=[G.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -G.b*y(3)+y(1)*y(2)-G.b*G.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = Modify_pointset(D,G) 

    tic;
    D.keys = keys(D.P);
    D.m = numEntries(D.P);
    disp("Modify Pointset - Initializing Variables: " + string(toc));

    key_time = 0; state_time = 0; iskey_time = 0;
    for l=1:D.m    % Check/Create Neighbors of Big Cells
        if(D.P(D.keys(l))>G.thresh)
            tic;
            current_key = D.keys(l); og_state = current_key{1}; key_time = key_time + toc;
    
            %   Faces (6)
            tic;
            new_key = {[og_state(1) + 1 og_state(2) og_state(3)]};
            state_time = state_time + toc;
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) og_state(3)]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) + 1 og_state(3)]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) - 1 og_state(3)]}; 
            state_time = state_time + toc;
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end  
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) og_state(3) + 1]}; 
            state_time = state_time + toc;
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end  
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) og_state(3) - 1]}; 
            state_time = state_time + toc;
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
    
            %   Edges (12)
            tic;
            new_key = {[og_state(1) + 1 og_state(2) + 1 og_state(3)]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) + 1 og_state(2) - 1 og_state(3)]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end  
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) + 1 og_state(3)]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) - 1 og_state(3)]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) + 1 og_state(2) og_state(3) + 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) + 1 og_state(2) og_state(3) - 1]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) og_state(3) + 1]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) og_state(3) - 1]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) + 1 og_state(3) + 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) + 1 og_state(3) - 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) - 1 og_state(3) + 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) og_state(2) - 1 og_state(3) - 1]}; state_time = state_time + toc; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;

            %   Corners (8)
            tic;
            new_key = {[og_state(1) + 1 og_state(2) + 1 og_state(3) + 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) + 1 og_state(2) + 1 og_state(3) - 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic; 
            new_key = {[og_state(1) + 1 og_state(2) - 1 og_state(3) + 1]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) + 1 og_state(3) + 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            tic;
            new_key = {[og_state(1) + 1 og_state(2) - 1 og_state(3) - 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) - 1 og_state(3) + 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) + 1 og_state(3) - 1]};
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
            tic;
            new_key = {[og_state(1) - 1 og_state(2) - 1 og_state(3) - 1]}; 
            state_time = state_time + toc; 
            tic;
            if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            iskey_time = iskey_time + toc;
        end
    end
    disp("Modify Pointset - Check/Create Neighbors - Key time: " + string(key_time));
    disp("Modify Pointset - Check/Create Neighbors - State converison: " + string(state_time));
    disp("Modify Pointset - Check/Create Neighbors - isKey Time: " + string(iskey_time));
    tic;
    D.keys = keys(D.P); D.values = values(D.P);
    D.n = numEntries(D.P);
    disp("Modify Pointset: Resetting values: " + string(toc));
    D = Update_vuw(D,G,D.m+1); 

    tic;
    D.keys = keys(D.P); D.values = values(D.P);
    v1 = D.v{1}; v2 = D.v{2}; v3 = D.v{3}; 
    u1 = D.u{1}; u2 = D.u{2}; u3 = D.u{3}; 
    w1 = D.w{1}; w2 = D.w{2}; w3 = D.w{3}; 
    for l=1:D.m                     % Remove Small Elements
        if(D.values(l) < G.thresh)
            D.P(D.keys(l)) = []; 
            v1(D.keys(l)) = []; v2(D.keys(l)) = []; v3(D.keys(l)) = [];
            u1(D.keys(l)) = []; u2(D.keys(l)) = []; u3(D.keys(l)) = [];
            w1(D.keys(l)) = []; w2(D.keys(l)) = []; w3(D.keys(l)) = [];
        end
    end
    disp("Modify Pointset - Removing empties: " + string(toc));
    tic;
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P);
    D.v = {v1, v2, v3}; D.u = {u1 u2 u3}; D.w = {w1, w2, w3};
     
    D.P(D.keys(1)) = 0;
    D.P(D.keys(:)) = max(D.P(D.keys(:)), 0);
    prob_sum = sum(values(D.P));
    D.P(D.keys(:)) = D.P(D.keys(:))./prob_sum;
    disp("Modify Pointset - Update Probabilities: " + string(toc));
end                                                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rotate_Plot(D,G,ys)                      
    N=round(2*G.L/G.dx)+1; M=(N-1)/2+1; Pfull=zeros(N,N,N); D.n = numEntries(D.P); D.keys = keys(D.P);
    for l=1:D.n, current_key = D.keys(l); state = current_key{1}; i = state(1)+M; j = state(2)+M; k = state(3)+M;
        if i>0 && i<=N && j>0 && j<=N && k>0 && k<=N, Pfull(j,i,k)=D.values(l); end, end
    figure(1)
    isosurface(-G.L:G.dx:G.L,-G.L:G.dx:G.L,-G.L:G.dx:G.L,Pfull,0.005); 
    isosurface(-G.L:G.dx:G.L,-G.L:G.dx:G.L,-G.L:G.dx:G.L,Pfull,0.0007); 
    isosurface(-G.L:G.dx:G.L,-G.L:G.dx:G.L,-G.L:G.dx:G.L,Pfull,0.0001); alpha(.5),
    colormap(cool); axis([-G.L G.L -G.L G.L -G.L G.L]);
    plot3(ys(1,:),ys(2,:),ys(3,:),'k-','linewidth',2);
    plot3(ys(1,end),ys(2,end),ys(3,end),'k*','linewidth',2); 
    axis equal; axis([-15 15 -25 25 -30 20]); drawnow;
    
    figure(2)
    isosurface(-G.L:G.dx:G.L,-G.L:G.dx:G.L,-G.L:G.dx:G.L,Pfull,0.0001); alpha(.5),
    colormap(cool); axis([-G.L G.L -G.L G.L -G.L G.L]);
    plot3(ys(1,:),ys(2,:),ys(3,:),'k-','linewidth',2);
    plot3(ys(1,end),ys(2,end),ys(3,end),'k*','linewidth',2);
    axis equal; axis([-15 15 -25 25 -30 20]); drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=RHS_P(D,G)
    tic;
    D.keys = keys(D.P); D.n = numEntries(D.P); D.P(D.keys(1)) = 0; K(1:D.n,1)=0; 
    disp("RHS_P - Initialize Variables: " + string(toc));
    % -------------------------------- NONCONSERVATIVE FORM -------------------------------
    % for l=2:D.n
    %   for d=1:G.d
    %       K(l,1)=K(l,1)-D.w(D.i(l,d),d)*(P(l)-P(D.i(l,d)))/G.dx-D.u(l,d)*(P(D.k(l,d))-P(l))/G.dx;    
    %   end 
    % end
    % --------------------------------- CONSERVATIVE FORM ---------------------------------    
    tic;
    for d=1:G.d 
        u = D.u{d}; w = D.w{d}; f = D.f{d};     
        for l=1:D.n
            current_key = D.keys(l); current_state = current_key{1};
            next_state = current_state; next_state(d) = next_state(d) + 1;
            next_key = {[next_state(1) next_state(2) next_state(3)]};
            if (isKey(D.P, next_key))
                f(current_key) = w(current_key)*D.P(current_key) + u(current_key)*D.P(next_key);
            else
                f(current_key) = w(current_key)*D.P(current_key);
            end
        end
        D.f{d} = f;
    end
    disp("RHS_P - Initialize Flux: " + string(toc));
    tic;
    for d=1:G.d
        vd = D.v{d}; ud = D.u{d}; wd = D.w{d}; fd = D.f{d};
        for l=1:D.n
            current_key = D.keys(l); current_state = current_key{1};
            prev_state = current_state; prev_state(d) = prev_state(d) - 1; prev_key = {[prev_state(1) prev_state(2) prev_state(3)]}; 
            if(isKey(D.P, prev_key))&&((D.P(current_key)>G.thresh)||(D.P(prev_key)>G.thresh))
                F=G.dt*(D.P(current_key)-D.P(prev_key))/(2*G.dx); 
                for e=1:G.d
                    if e~=d
                        ue = D.u{e}; we = D.w{e}; fe=D.f{e};
                        fe(current_key) = fe(current_key) - we(current_key)*wd(prev_key)*F; 
                        fe(prev_key) = fe(prev_key) - we(prev_key)*ud(prev_key)*F;
                        new_prev_state = current_state; new_prev_state(e) = new_prev_state(e) - 1; new_prev_key = {[new_prev_state(1) new_prev_state(2) new_prev_state(3)]};
                        if(isKey(D.P, new_prev_key))
                            fe(new_prev_key) = fe(new_prev_key) - ue(new_prev_key)*wd(prev_key)*F;
                        end
                        new_prev_state = prev_state; new_prev_state(e) = new_prev_state(e) - 1; new_prev_key = {[new_prev_state(1) new_prev_state(2) new_prev_state(3)]};
                        if(isKey(D.P, new_prev_key))
                            fe(new_prev_key) = fe(new_prev_key) - ue(new_prev_key)*ud(prev_key)*F;
                        end
                        D.u{e} = ue; D.w{e} = we; D.f{e} = fe;
                    end
                end
                if vd(prev_key)>0
                    new_prev_state = prev_state; new_prev_state(d) = new_prev_state(d) - 1;
                    new_prev_key = {[new_prev_state(1) new_prev_state(2) new_prev_state(3)]};
                    if(isKey(D.P, new_prev_key))
                        th = (D.P(prev_key) - D.P(new_prev_key))/(D.P(current_key)-D.P(prev_key));
                    else
                        th = D.P(prev_key)/(D.P(current_key)-D.P(prev_key));
                    end
                else
                    next_state = current_state; next_state(d) = next_state(d) + 1;
                    next_key = {[next_state(1) next_state(2) next_state(3)]};
                    if(isKey(D.P, next_key))
                        th = (D.P(next_key) - D.P(current_key))/(D.P(current_key) - D.P(prev_key));
                    else
                        th = -D.P(current_key)/(D.P(current_key) - D.P(prev_key));
                    end
                end
                t=abs(vd(prev_key)); fd(prev_key) = fd(prev_key) + t*(G.dx/G.dt-t)*F*MC(th);% Flux: use MC or VL.
            end
        end
    end
    disp("RHS_P - Solve v,u,w,th,t: " + string(toc));
    tic;
    for l=1:D.n
        for d=1:G.d
            fd = D.f{d}; current_key = D.keys(l); current_state = current_key{1};
            prev_state = current_state; prev_state(d) = prev_state(d) - 1;
            prev_key = {[prev_state(1) prev_state(2) prev_state(3)]};
            if(isKey(D.P, prev_key))
                K(l,1)=K(l,1)-(fd(current_key)-fd(prev_key))/G.dx; 
            end
        end
    end
    disp("RHS_P - Solve K: " + string(toc));
end % function RHS_P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi]=MC(th), phi=max(0,min([(1+th)/2 2 2*th]));   end   % Flux limiters
function [phi]=VL(th), phi=min((th+abs(th))/(1+abs(th)),0); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function match = keys_match(dict_list)
    match = 1;
    for i=1:length(dict_list)
        for j=1:length(dict_list)
            if i~=j
                match = match*(isequal(keys(dict_list{i}), keys(dict_list{j})));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = isEdge(D, G, state, d0)
    bool = 1;
    for d=1:G.d
        prev_state = state; prev_state(d) = prev_state(d) - 1;
        bool = bool*(isKey(D.P, state_conversion(prev_state)));
        if d==d0
            prev_state(d) = prev_state(d) - 1;
            bool = bool*(isKey(D.P, state_conversion(prev_state)));
        end
    end
    bool = ~bool;
end