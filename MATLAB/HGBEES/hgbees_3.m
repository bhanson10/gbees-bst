hD = hgbees;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hD = hgbees
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
T=1; G.thresh=0.00002; G.start=[-11.5; -10; 9.5]; G.max = 10000;
G.dt=.0005; dt=.005; G.dx=0.4; G.d=3; G.sigma=4; G.b=1; G.r=48; G.L=30;
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%% 
G.Y = eye(G.d,'int16'); [hD] = Initialize_D(G);

y=G.start; ys=y; figure(1); clf;
for timestep=1:10000
  k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
  ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
end
plot3(ys(1,:),ys(2,:),ys(3,:),'g-','linewidth',1); view(-109,14);  hold on;
lighting phong; light('Position',[-1 0 0]); drawnow;

G.xrange = [min(ys(1,:))-1, max(ys(1,:))+1]; 
G.yrange = [min(ys(2,:))-1, max(ys(2,:))+1];
G.zrange = [min(ys(3,:))-1, max(ys(3,:))+1];
[hD] = Initialize_Neighbors(G);

figure(2); clf; y=G.start; ys=y; 
for timestep=1:T/dt
  k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
  ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
end
plot3(ys(1,:),ys(2,:),ys(3,:),'k-','linewidth',2); view(-109,14);  hold on;
plot3(ys(1,1),ys(2,1),ys(3,1),'k*'), plot3(ys(1,end),ys(2,end),ys(3,end),'k*');
lighting phong; light('Position',[-1 0 0]); drawnow;

for P=1:200; y=G.start+0.5*randn(3,1); ys=y;
  for timestep=1:T/dt
    k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
    ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
  end
  plot3(ys(1,:),ys(2,:),ys(3,:),'c-.','linewidth',0.3); 
  plot3(ys(1,1),ys(2,1),ys(3,1),'k+'), plot3(ys(1,end),ys(2,end),ys(3,end),'k+');
  plot3(ys(1,41),ys(2,41),ys(3,41),'k+'); plot3(ys(1,81),ys(2,81),ys(3,81),'k+');
  plot3(ys(1,121),ys(2,121),ys(3,121),'k+'); plot3(ys(1,161),ys(2,161),ys(3,161),'k+');
  drawnow;
end


y=G.start; ys=y; t=0; [hD]=Modify_pointset(hD,G); Rotate_Plot(hD,G,ys); pause(inf);
%for timestep=1:T/G.dt
for timestep=1:10, disp(timestep); t=t+G.dt; if mod(timestep,1)==0, [hD]=Modify_pointset(hD,G); end
  K=RHS_P(hD,G); hD.values = values(hD.P); hD.keys = keys(hD.P); hD.P(hD.keys) = hD.values + G.dt.*K;
  k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
  ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew;
  %if mod(timestep,400)==0, Rotate_Plot(hD,G,ys),  end
end, %Rotate_Plot(hD,G,ys),

%{
figure(1);
view(-109,14); print -depsc2 -image -r600 pdfA.v1.eps
view(-31,2);   print -depsc2 -image -r600 pdfA.v2.eps

figure(2);
view(-109,14); print -depsc2 -image -r600 trajA.v1.eps
view(-31,2);   print -depsc2 -image -r600 trajA.v2.eps
%}

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = Initialize_D(G)
    D.P = dictionary();
    D.v=zeros(G.max,G.d); D.w=D.v; D.u=D.v; D.f=D.v; 
   
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
    D.m=numEntries(D.P); D.n=D.m; D = Update_vuw(D,G,1); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=Update_vuw(D,G,b)
    for l=b:D.n
      current_key = D.keys(l);
      x=G.dx.*(key_conversion(current_key)); xh=G.dx/2;
      D.v(l,1)=G.sigma*(x(2)-(x(1)+xh));
      D.v(l,2)=-(x(2)+xh)-x(1)*x(3); 
      D.v(l,3)=-G.b*(x(3)+xh)+x(1)*x(2)-G.b*G.r;
      % ------------------------- (keep one of the above 2 sections) ----------------------
      for d=1:G.d, D.w(l,d)=max(D.v(l,d),0); D.u(l,d)=min(D.v(l,d),0); end  % Init u and w.
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
    D.eys = keys(D.P); D.m = numEntries(D.P);
    for l=1:D.m    % Check/Create Neighbors of Big Cells
        if(D.P(D.keys(l))>G.thresh)
            key_neighbors = D.neighbors(D.keys(l));
            for i=1:G.d
                slice = key_neighbors{i};
                for j=1:G.d
                    for k=1:G.d
                        new_key = slice(j,k);
                        if(~isKey(D.P, new_key)), D.P(new_key) = 0;end
                    end
                end
            end               
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P);
    D = Update_vuw(D,G,D.m+1); 
    for l=1:D.m                     % Remove Small Elements
        if(D.values(l) < G.thresh)
            D.P(D.keys(l)) = []; D.v(D.keys(l)) = []; D.u(D.keys(l)); D.w(D.keys(l)) = [];
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P);    
    D.P(D.keys(1)) = 0; D.P(D.keys(:)) = max(D.P(D.keys(:)), 0);
    prob_sum = sum(values(D.P)); D.P(D.keys(:)) = D.P(D.keys(:))./prob_sum;
end                                                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rotate_Plot(D,G,ys)                      
    N=round(2*G.L/G.dx)+1; M=(N-1)/2+1; Pfull=zeros(N,N,N); D.n = numEntries(D.P); D.keys = keys(D.P);
    for l=1:D.n, state = key_conversion(D.keys(l)); i = state(1)+M; j = state(2)+M; k = state(3)+M;
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
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P); D.P(D.keys(1)) = 0; K(1:D.n,1)=0; 

    for l=1:D.n
        current_key = D.keys(l); current_state = key_conversion(current_key);
        next_state = current_state; next_state(d) = next_state(d) + 1;
        next_key = state_conversion(next_state);
        current_P = D.P(current_key); next_p = D.p(next_key);
        u = D.u(current_key); u = u{1}; w = D.w(current_key); w = w{1};
        if (isKey(D.P, next_key))
            D.f(current_key) = {[w(1)*current_P + u(1)*next_P w(2)*current_P + u(2)*next_P w(3)*current_P + u(3)*next_P]};
        else
            D.f(current_key) = {[w(1)*current_P w(2)*current_P w(3)*current_P]};
        end
    end
    for d=1:G.d
        vd = D.v{d}; ud = D.u{d}; wd = D.w{d}; fd = D.f{d};
        for l=1:D.n
            current_key = D.keys(l); current_state = key_conversion(current_key);
            prev_state = current_state; prev_state(d) = prev_state(d) - 1; prev_key = state_conversion(prev_state); 
            if(isKey(D.P, prev_key))&&((D.P(current_key)>G.thresh)||(D.P(prev_key)>G.thresh))
                F=G.dt*(D.P(current_key)-D.P(prev_key))/(2*G.dx); 
                for e=1:G.d
                    if e~=d
                        ue = D.u{e}; we = D.w{e}; fe=D.f{e};
                        fe(current_key) = fe(current_key) - we(current_key)*wd(prev_key)*F; 
                        fe(prev_key) = fe(prev_key) - we(prev_key)*ud(prev_key)*F;
                        new_prev_state = current_state; new_prev_state(e) = new_prev_state(e) - 1; new_prev_key = state_conversion(new_prev_state);
                        if(isKey(D.P, new_prev_key))
                            fe(new_prev_key) = fe(new_prev_key) - ue(new_prev_key)*wd(prev_key)*F;
                        end
                        new_prev_state = prev_state; new_prev_state(e) = new_prev_state(e) - 1; new_prev_key = state_conversion(new_prev_state);
                        if(isKey(D.P, new_prev_key))
                            fe(new_prev_key) = fe(new_prev_key) - ue(new_prev_key)*ud(prev_key)*F;
                        end
                        D.u{e} = ue; D.w{e} = we; D.f{e} = fe;
                    end
                end
                if vd(prev_key)>0
                    new_prev_state = prev_state; new_prev_state(d) = new_prev_state(d) - 1;
                    new_prev_key = state_conversion(new_prev_state);
                    if(isKey(D.P, new_prev_key))
                        th = (D.P(prev_key) - D.P(new_prev_key))/(D.P(current_key)-D.P(prev_key));
                    else
                        th = D.P(prev_key)/(D.P(current_key)-D.P(prev_key));
                    end
                else
                    next_state = current_state; next_state(d) = next_state(d) + 1;
                    next_key = state_conversion(next_state);
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
    
    for l=1:D.n
        for d=1:G.d
            fd = D.f{d}; current_key = D.keys(l); current_state = key_conversion(current_key);
            prev_state = current_state; prev_state(d) = prev_state(d) - 1;
            prev_key = state_conversion(prev_state);
            if(isKey(D.P, prev_key))
                K(l,1)=K(l,1)-(fd(current_key)-fd(prev_key))/G.dx; 
            end
        end
    end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = Initialize_Neighbors(G) 
    all_keys = {}; x_count = 0; y_count = 0;
    for k=round((G.zrange(1))/G.dx):round((G.zrange(2))/G.dx)
        for i=round((G.xrange(1))/G.dx):round((G.xrange(2))/G.dx)
            x_count = x_count + 1;
            for j=round((G.yrange(1))/G.dx):round((G.yrange(2))/G.dx)
                y_count = y_count + 1;
                state = [i j k]; key = state_conversion(state);
                matrix(x_count, y_count) = string(key);
            end
            y_count = 0;
        end
        x_count = 0;
        all_keys{end+1} = matrix;
    end

    D.neighbors = containers.Map();
    z_width = (round((G.zrange(2))/G.dx)-round((G.zrange(1))/G.dx))+1;
    y_width = (round((G.yrange(2))/G.dx)-round((G.yrange(1))/G.dx))+1;
    x_width = (round((G.xrange(2))/G.dx)-round((G.xrange(1))/G.dx))+1;
    for i=1:z_width
        slice = all_keys{i};
        for k=1:x_width
            for j=1:y_width
                n_key = slice(k,j); 
    
                %Back Z
                if(i==1)
                    back_slice = ["0" "0" "0"; "0" "0" "0"; "0" "0" "0"];
                else
                    slice = all_keys{i-1};
                    if(k==1)&&(j==1)
                        back_slice = ["0" "0" "0"; "0" slice(k,j) slice(k,j+1); "0" slice(k+1,j) slice(k+1,j+1)];                
                    elseif(k==x_width)&&(j==y_width)
                        back_slice = [slice(k-1, j-1) slice(k-1, j) "0"; slice(k, j-1) slice(k, j) "0"; "0" "0" "0"];                
                    elseif(k==1)&&(j==y_width)
                        back_slice = ["0" "0" "0"; slice(k, j-1) slice(k,j) "0"; slice(k+1,j-1) slice(k+1, j) "0"];
                    elseif(k==x_width)&&(j==1)
                        back_slice = ["0" slice(k-1, j) slice(k-1, j+1); "0" slice(k,j) slice(k, j+1); "0" "0" "0"];
                    elseif(k==1)
                        back_slice = ["0" "0" "0"; slice(k, j-1) slice(k,j) slice(k,j+1); slice(k+1, j-1) slice(k+1,j) slice(k+1,j+1)];
                    elseif(j==1)
                        back_slice = ["0" slice(k-1,j) slice(k-1,j+1); "0" slice(k,j) slice(k,j+1); "0" slice(k+1,j) slice(k+1,j+1)];                   
                    elseif(k==x_width)
                        back_slice = [slice(k-1,j-1) slice(k-1, j) slice(k-1, j+1); slice(k,j-1) slice(k, j) slice(k, j+1); "0" "0" "0"];
                    elseif(j==y_width)
                        back_slice = [slice(k-1,j-1) slice(k-1,j) "0"; slice(k,j-1) slice(k,j) "0"; slice(k+1,j-1) slice(k+1,j) "0"];
                    else
                        back_slice = [slice(k-1,j-1) slice(k-1,j) slice(k-1,j+1); slice(k,j-1) slice(k,j) slice(k,j+1); slice(k+1,j-1) slice(k+1,j) slice(k+1,j+1)];
                    end
                end
    
                %Current z
                slice = all_keys{i};
                if(k==1)&&(j==1)
                    cur_slice = ["0" "0" "0"; "0" slice(k,j) slice(k,j+1); "0" slice(k+1,j) slice(k+1,j+1)];                
                elseif(k==x_width)&&(j==y_width)
                    cur_slice = [slice(k-1, j-1) slice(k-1, j) "0"; slice(k, j-1) slice(k, j) "0"; "0" "0" "0"];                
                elseif(k==1)&&(j==y_width)
                    cur_slice = ["0" "0" "0"; slice(k, j-1) slice(k,j) "0"; slice(k+1,j-1) slice(k+1, j) "0"];
                elseif(k==x_width)&&(j==1)
                    cur_slice = ["0" slice(k-1, j) slice(k-1, j+1); "0" slice(k,j) slice(k, j+1); "0" "0" "0"];
                elseif(k==1)
                    cur_slice = ["0" "0" "0"; slice(k, j-1) slice(k,j) slice(k,j+1); slice(k+1, j-1) slice(k+1,j) slice(k+1,j+1)];
                elseif(j==1)
                    cur_slice = ["0" slice(k-1,j) slice(k-1,j+1); "0" slice(k,j) slice(k,j+1); "0" slice(k+1,j) slice(k+1,j+1)];                   
                elseif(k==x_width)
                    cur_slice = [slice(k-1,j-1) slice(k-1, j) slice(k-1, j+1); slice(k,j-1) slice(k, j) slice(k, j+1); "0" "0" "0"];
                elseif(j==y_width)
                    cur_slice = [slice(k-1,j-1) slice(k-1,j) "0"; slice(k,j-1) slice(k,j) "0"; slice(k+1,j-1) slice(k+1,j) "0"];
                else
                    cur_slice = [slice(k-1,j-1) slice(k-1,j) slice(k-1,j+1); slice(k,j-1) slice(k,j) slice(k,j+1); slice(k+1,j-1) slice(k+1,j) slice(k+1,j+1)];
                end
    
                %Forward Z
                if(i==z_width)
                    for_slice = ["0" "0" "0"; "0" "0" "0"; "0" "0" "0"];
                else
                    slice = all_keys{i+1};
                    if(k==1)&&(j==1)
                        for_slice = ["0" "0" "0"; "0" slice(k,j) slice(k,j+1); "0" slice(k+1,j) slice(k+1,j+1)];                
                    elseif(k==x_width)&&(j==y_width)
                        for_slice = [slice(k-1, j-1) slice(k-1, j) "0"; slice(k, j-1) slice(k, j) "0"; "0" "0" "0"];                
                    elseif(k==1)&&(j==y_width)
                        for_slice = ["0" "0" "0"; slice(k, j-1) slice(k,j) "0"; slice(k+1,j-1) slice(k+1, j) "0"];
                    elseif(k==x_width)&&(j==1)
                        for_slice = ["0" slice(k-1, j) slice(k-1, j+1); "0" slice(k,j) slice(k, j+1); "0" "0" "0"];
                    elseif(k==1)
                        for_slice = ["0" "0" "0"; slice(k, j-1) slice(k,j) slice(k,j+1); slice(k+1, j-1) slice(k+1,j) slice(k+1,j+1)];
                    elseif(j==1)
                        for_slice = ["0" slice(k-1,j) slice(k-1,j+1); "0" slice(k,j) slice(k,j+1); "0" slice(k+1,j) slice(k+1,j+1)];                   
                    elseif(k==x_width)
                        for_slice = [slice(k-1,j-1) slice(k-1, j) slice(k-1, j+1); slice(k,j-1) slice(k, j) slice(k, j+1); "0" "0" "0"];
                    elseif(j==y_width)
                        for_slice = [slice(k-1,j-1) slice(k-1,j) '0'; slice(k,j-1) slice(k,j) '0'; slice(k+1,j-1) slice(k+1,j) '0'];
                    else
                        for_slice = [slice(k-1,j-1) slice(k-1,j) slice(k-1,j+1); slice(k,j-1) slice(k,j) slice(k,j+1); slice(k+1,j-1) slice(k+1,j) slice(k+1,j+1)];
                    end
                end
                D.neighbors(n_key) = {back_slice, cur_slice, for_slice};
            end
        end
    end
end