function hgbees_bewley
%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%% 
T=1; G.thresh=0.00002; G.start=[-11.5; -10; 9.5]; G.dt=.0005;
dt=.005; G.dx=0.4; G.d=3; G.sigma=4; G.b=1; G.r=48; G.L=30; G.xh=G.dx/2;
G.N_bits = 8; G.N_data = G.d; fac = uint64(2^G.N_bits); for i=1:G.N_data, G.fac(i) = fac^(i-1); end 
G.offset32=int32(fac/2); G.offset64=int64(G.offset32);
%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%% 
[hD] = Initialize_D(G); 

figure(1); clf; y=G.start; ys=y; 
for timestep=1:10000
  k1=RHS(y,G); k2=RHS(y+(dt/2)*k1,G); k3=RHS(y+(dt/2)*k2,G); k4=RHS(y+dt*k3,G);    
  ynew=y+(dt/6)*k1+(dt/3)*(k2+k3)+(dt/6)*k4; ys=[ys ynew]; y=ynew;
end
plot3(ys(1,:),ys(2,:),ys(3,:),'g-','linewidth',1); view(-109,14);  hold on;
lighting phong; light('Position',[-1 0 0]); drawnow;

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

y=G.start; ys=y; t=0; [hD]=Modify_pointset(hD,G); %Rotate_Plot(hD,G,ys); 

%for timestep=1:T/G.dt
for timestep=1:T/G.dt, disp("Timestep: " + string(timestep)); t=t+G.dt; if mod(timestep,1)==0, [hD]=Modify_pointset(hD,G); end
  K=RHS_P(hD,G); hD.keys = keys(hD.P); hD.P(hD.keys) = hD.P(hD.keys) + G.dt.*K;
  k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
  ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew; 
  if mod(timestep,25)==0, Rotate_Plot(hD,G,ys),  end
end, Rotate_Plot(hD,G,ys),

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
    D.P = dictionary(); D.P(uint64(0)) = 0;
    D.v = dictionary(); D.v(uint64(0)) = {zeros(1,G.d)};
    D.u = D.v; D.w = D.v; D.f = D.v;
   
    for i=round((G.start(1)-2)/G.dx):round((G.start(1)+2)/G.dx)
        for j=round((G.start(2)-2)/G.dx):round((G.start(2)+2)/G.dx)
            for k=round((G.start(3)-2)/G.dx):round((G.start(3)+2)/G.dx)
                state = [i j k]; key = state_conversion(state,G);
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
        current_key = D.keys(l); current_state = double(key_conversion(current_key,G));
        x=G.dx.*current_state; 
        v1=G.sigma*(x(2)-(x(1)+G.xh)); v2=-(x(2)+G.xh)-x(1)*x(3); v3=-G.b*(x(3)+G.xh)+x(1)*x(2)-G.b*G.r; 
        D.v(current_key) = {[v1 v2 v3]};
        D.u(current_key)={[min(v1, 0) min(v2, 0) min(v3, 0)]};
        D.w(current_key)={[max(v1, 0) max(v2, 0) max(v3, 0)]};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(y,G)                          
    f=[G.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -G.b*y(3)+y(1)*y(2)-G.b*G.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = state_conversion(state,G)
    state = int32(state); state = uint64(state + G.offset32);
    key=uint64(0);
    for i=1:G.N_data; key=key+state(i)*G.fac(i); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = key_conversion(key,G)
    state = zeros(1,G.d);
    for i=G.N_data:-1:1
        state(i)=idivide(key,G.fac(i),'floor');
        key=key-state(i)*G.fac(i);
    end
    state=int64(state)-G.offset64;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = Modify_pointset(D,G) 
    D.keys = keys(D.P); D.m = numEntries(D.P);
    
    for l=2:D.m                       % Check/Create Neighbors of Big Cells
        if(D.P(D.keys(l))>=G.thresh)
            current_key = D.keys(l); state = key_conversion(current_key, G);
            x_neighbors = state(1)-1:state(1)+1; 
            y_neighbors = state(2)-1:state(2)+1;
            z_neighbors = state(3)-1:state(3)+1;

            [X,Y,Z] = meshgrid(x_neighbors,y_neighbors,z_neighbors); 
            X = reshape(X,[1,3^G.d]); Y = reshape(Y,[1,3^G.d]); Z = reshape(Z,[1,3^G.d]);
        
            for i=1:3^G.d
                new_state = [X(i) Y(i) Z(i)]; new_key = state_conversion(new_state, G);
                if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            end   
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P);
    D = Initialize_vuw(D,G,D.m+1); 
    for l=2:D.m                                     % Remove Small Elements
        if(D.values(l) < G.thresh)&&(no_neighbors(D,G,l))
            D.P(D.keys(l)) = []; D.v(D.keys(l)) = []; D.u(D.keys(l)) = []; D.w(D.keys(l)) = [];
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P);    
    D.P(D.keys) = max(D.P(D.keys), 0); D.values = values(D.P); 
    prob_sum = sum(values(D.P)); D.P(D.keys) = D.values./prob_sum; D.values = values(D.P);
end                                                                                                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rotate_Plot(D,G,ys)                      
    N=round(2*G.L/G.dx)+1; M=(N-1)/2+1; Pfull=zeros(N,N,N); D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P);
    for l=2:D.n, state = key_conversion(D.keys(l),G); i = state(1)+M; j = state(2)+M; k = state(3)+M;
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
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P); K(1:D.n,1)=0; 
    D.f(D.keys) = {zeros(1,G.d)}; states = {};
    for l=2:D.n
        l_key = D.keys(l); l_state = key_conversion(l_key,G); states{end+1} = l_state;
        f_l = D.f(l_key); f_l = f_l{1}; u_l = D.u(l_key); u_l = u_l{1}; w_l = D.w(l_key); w_l = w_l{1};  
        for d=1:G.d
            k_state = l_state; k_state(d) = k_state(d)+1; k_key = state_conversion(k_state,G);
            if(~isKey(D.P,k_key)),k_key = uint64(0);end
            f_l(d) = w_l(d)*D.P(l_key) + u_l(d)*D.P(k_key);
        end
        D.f(l_key) = {f_l};
    end

    for d=1:G.d
        for l=2:D.n
            l_key = D.keys(l); l_state = states{l-1};
            f_l = D.f(l_key); f_l = f_l{1}; w_l = D.w(l_key); w_l = w_l{1}; 
            i_state = l_state; i_state(d) = i_state(d)-1; i_key = state_conversion(i_state,G);
            if(~isKey(D.P,i_key)),i_key = uint64(0);end
            f_i = D.f(i_key); f_i = f_i{1}; w_i = D.w(i_key); w_i = w_i{1}; u_i = D.u(i_key); u_i = u_i{1};
            v_i = D.v(i_key); v_i = v_i{1};
            if (D.P(l_key)>G.thresh)||(D.P(i_key)>G.thresh)
                F=G.dt*(D.P(l_key)-D.P(i_key))/(2*G.dx);
                for e=1:G.d
                    if e~=d
                        j_state = l_state; j_state(e) = j_state(e)-1; j_key = state_conversion(j_state,G);
                        if(~isKey(D.P,j_key)),j_key = uint64(0);end
                        f_j = D.f(j_key); f_j = f_j{1}; w_j = D.w(j_key); w_j = w_j{1}; u_j = D.u(j_key); u_j = u_j{1};
                        p_state = i_state; p_state(e) = p_state(e)-1; p_key = state_conversion(p_state,G);
                        if(~isKey(D.P,p_key)),p_key = uint64(0);end
                        f_p = D.f(p_key); f_p = f_p{1}; w_p = D.w(p_key); w_p = w_p{1}; u_p = D.u(p_key); u_p = u_p{1};

                        f_l(e) = f_l(e)-w_l(e)*w_i(d)*F;
                        f_j(e) = f_j(e)-u_j(e)*w_i(d)*F; D.f(j_key) = {f_j};
                        f_i(e) = f_i(e)-w_i(e)*u_i(d)*F; 
                        f_p(e) = f_p(e)-u_p(e)*u_i(d)*F; D.f(p_key) = {f_p};
                    end
                end
                D.f(l_key) = {f_l};
                D.f(i_key) = {f_i};
                
                i_i_state = i_state; i_i_state(d) = i_i_state(d)-1; i_i_key = state_conversion(i_i_state,G);
                if(~isKey(D.P,i_i_key)),i_i_key = uint64(0);end
                k_state = l_state; k_state(d) = k_state(d)+1; k_key = state_conversion(k_state,G);
                if(~isKey(D.P,k_key)),k_key = uint64(0);end
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
    
    for l=2:D.n
        l_key = D.keys(l); l_state = states{l-1};
        f_l = D.f(l_key); f_l = f_l{1};
        for d=1:G.d
            i_state = l_state; i_state(d) = i_state(d)-1; i_key = state_conversion(i_state,G);
            if(~isKey(D.P,i_key)),i_key = uint64(0);end
            f_i = D.f(i_key); f_i = f_i{1};
            K(l,1)=K(l,1)-(f_l(d)-f_i(d))/G.dx;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi]=MC(th), phi=max(0,min([(1+th)/2 2 2*th]));   end   % Flux limiters
function [phi]=VL(th), phi=min((th+abs(th))/(1+abs(th)),0); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool=no_neighbors(D,G,l)
    bool = 1;
    current_key = D.keys(l); current_state = key_conversion(current_key,G);

    x_neighbors = current_state(1)-1:current_state(1)+1;
    y_neighbors = current_state(2)-1:current_state(2)+1;
    z_neighbors = current_state(3)-1:current_state(3)+1;
    
    [X,Y,Z] = ndgrid(x_neighbors,y_neighbors,z_neighbors); 
    X = reshape(X,[1,3^G.d]); Y = reshape(Y,[1,3^G.d]); Z = reshape(Z,[1,3^G.d]);

    for i=1:3^G.d
        new_state = [X(i) Y(i) Z(i)]; new_key = state_conversion(new_state, G);
        if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%