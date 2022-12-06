% Used for testing how long each GBEES method takes

%GBEES
disp("%%%%%%%%%%% ENTERING GBEES %%%%%%%%%%%");
disp(" ");

T=1; G.thresh=0.00002; G.max=10000; G.start=[-11.5; -10; 9.5]; 
G.dt=.0005; dt=.005; G.dx=0.4; G.d=3; G.sigma=4; G.b=1; G.r=48; G.L=30;

G.Y=eye(G.d,'int16'); [D]=G_Initialize_D(G); h1=round(G.dx/G.dt); 

y=G.start; ys=y; t=0; [D]=G_Modify_pointset(D,G); G_ind_time = []; G_total_time = []; G_size = [];
%for timestep=1:T/G.dt,
for timestep=1:50, disp("Timestep: " + string(timestep)); tic; if mod(timestep,1)==0, [D]=G_Modify_pointset(D,G); end
    K=G_RHS_P(D,D.P(1:D.n),G); D.P(2:D.n,1)=D.P(2:D.n)+G.dt*K(2:D.n);                
    
    k1=RHS(y,G); k2=RHS(y+(G.dt/2)*k1,G); k3=RHS(y+(G.dt/2)*k2,G); k4=RHS(y+G.dt*k3,G);    
    ynew=y+(G.dt/6)*k1+(G.dt/3)*(k2+k3)+(G.dt/6)*k4; ys=[ys ynew]; y=ynew;
    
    ind_time = toc;
    if timestep == 1
        step_time = ind_time;
    else
        step_time = G_total_time(end) + ind_time;
    end
    
    G_ind_time = [G_ind_time ind_time];
    G_total_time = [G_total_time step_time];
    G_size = [G_size D.n-1];

end


%HGBEES_1
disp("%%%%%%%%%%% ENTERING HGBEES_1 %%%%%%%%%%%");
disp(" ");

G.N_bits = 8; G.N_data = G.d; G.fac=uint64(2^G.N_bits); G.offset32=int32(G.fac/2); G.offset64=int64(G.offset32);
[hD] = H1_Initialize_D(G); 
t=0; [hD]=H1_Modify_pointset(hD,G); H1_ind_time = []; H1_total_time = []; H1_size = [];

for timestep=1:50, disp("Timestep: " + string(timestep)); tic; t=t+G.dt; if mod(timestep,1)==0, [hD]=H1_Modify_pointset(hD,G); end
    K=H1_RHS_P(hD,G); hD.keys = keys(hD.P); hD.P(hD.keys) = hD.P(hD.keys) + G.dt.*K;
    
    ind_time = toc;
    if timestep == 1
        step_time = ind_time;
    else
        step_time = H1_total_time(end) + ind_time;
    end
    
    H1_ind_time = [H1_ind_time ind_time];
    H1_total_time = [H1_total_time step_time];
    H1_size = [H1_size length(hD.keys)-1];
end

%HGBEES_2
disp("%%%%%%%%%%% ENTERING HGBEES_2 %%%%%%%%%%%");
disp(" ");

[hD] = H2_Initialize_D(G); 
y=G.start; ys=y; t=0; [hD]=H2_Modify_pointset(hD,G); H2_ind_time = []; H2_total_time = []; H2_size = [];

for timestep=1:50, disp("Timestep: " + string(timestep)); tic; t=t+G.dt; if mod(timestep,1)==0, [hD]=H2_Modify_pointset(hD,G); end
    K=H2_RHS_P(hD,G); hD.keys = keys(hD.P); hD.P(hD.keys) = hD.P(hD.keys) + G.dt.*K;
    
    ind_time = toc;
    if timestep == 1
        step_time = ind_time;
    else
        step_time = H2_total_time(end) + ind_time;
    end
    
    H2_ind_time = [H2_ind_time ind_time];
    H2_total_time = [H2_total_time step_time];
    H2_size = [H2_size length(hD.keys)-1];
end

figure(1); clf; hold on
timestep = [1:50];
plot(timestep, G_total_time, 'k', 'LineWidth', 1, 'DisplayName', 'GBEES');
plot(timestep, H1_total_time, 'b', 'LineWidth', 1, 'DisplayName', 'HGBEES1');
plot(timestep, H2_total_time, 'r', 'LineWidth', 1, 'DisplayName', 'HGBEES2');
title('Total time of GBEES Variants', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "northwest";
lgd.FontSize = 10;
xlabel('Step \#', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('time (s)', 'Interpreter', 'Latex', 'FontSize', 10)

figure(2); clf; hold on
scatter(G_size, G_ind_time, 'k', 'filled', 'DisplayName', 'GBEES');
scatter(H1_size, H1_ind_time, 'b', 'filled', 'DisplayName', 'HGBEES1');
scatter(H2_size, H2_ind_time, 'r', 'filled', 'DisplayName', 'HGBEES2');
title('Time Complexity of GBEES Variants', 'Interpreter','Latex', 'FontSize', 14);
lgd = legend;
lgd.Location = "northwest";
lgd.FontSize = 10;
xlabel('Size of Dictionary', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('Time of step (s)', 'Interpreter', 'Latex', 'FontSize', 10)

disp("Total time of GBEES: " +string(G_total_time(end)));
disp("Total time of HBEES1: " +string(H1_total_time(end)));
disp("Total time of HGBEES2: " +string(H2_total_time(end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     GBEES FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=G_Initialize_D(G)
D.P=zeros(G.max,1);   D.j=zeros(G.max,G.d,'int16'); D.k=ones(G.max,G.d,'int16');
D.v=zeros(G.max,G.d); D.i=D.k; D.w=D.v; D.u=D.v; D.f=D.v; l=1;
% ------------------------------------- 3D Gaussian -----------------------------------
for i=round(-13.5/G.dx):round(-9.5/G.dx),
 for j=round(-12/G.dx):round(-8/G.dx),
  for k=round(7.5/G.dx):round(11.5/G.dx),
    x=(i*G.dx + 11.5)^2+(j*G.dx + 10)^2+(k*G.dx - 9.5)^2;
    l=l+1; D.j(l,1)=i; D.j(l,2)=j; D.j(l,3)=k; D.P(l)=exp(-4.*x/2.);
end, end, end
% ------------------------------------ Single Point -----------------------------------
% l=2; D.j(l,1:3)=round(9/G.dx); D.P(l)=1/G.dx^3;
% -------------------------------------------------------------------------------------
D.m=l; D.n=D.m; [D]=Initialize_vuwik(D,G,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=RHS(y,G)                          
f=[G.sigma*(y(2)-y(1));  -y(2)-y(1)*y(3);  -G.b*y(3)+y(1)*y(2)-G.b*G.r];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=Initialize_vuwik(D,G,b);
for l=b:D.n
  x(1:G.d)=G.dx*double(D.j(l,1:G.d)); xh=G.dx/2;
  D.v(l,1)=G.sigma*(x(2)-(x(1)+xh));
  D.v(l,2)=-(x(2)+xh)-x(1)*x(3); 
  D.v(l,3)=-G.b*(x(3)+xh)+x(1)*x(2)-G.b*G.r;
  % ------------------------- (keep one of the above 2 sections) ----------------------
  for d=1:G.d, D.w(l,d)=max(D.v(l,d),0); D.u(l,d)=min(D.v(l,d),0); end  % Init u and w.
end
D.i(b:D.n,:)=ones(D.n-b+1,G.d,'int16'); D.k(b:D.n,:)=ones(D.n-b+1,G.d,'int16'); 
for l=D.n:-1:b,                                         % Search list for neighbors to
  for t=2:l-1, diff=(D.j(l,:)-D.j(t,:));       % init i and k.  For large D.n,
    if sum(~diff)==G.d-1; [Y,d] = max(abs(diff));       % THIS IS THE EXPENSIVE BIT!
      if     D.j(l,d)==D.j(t,d)+1, D.i(l,d)=t; D.k(t,d)=l; % Due to careful programming
	  elseif D.j(t,d)==D.j(l,d)+1, D.k(l,d)=t; D.i(t,d)=l; % of Modify_pointset, we do
  end, end, end                                            % not call it very often...
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=G_Modify_pointset(D,G)               % Shift dataset to relevant gridpoints.
D.m=D.n; l=D.m;

while D.P(l)<G.thresh; l=l-1; D.m=D.m-1; end, l=l-1;  % Find last big element.
while l>1, if D.P(l)<G.thresh,     % We now move all big elements to the range 2:D.m,
  for d=1:G.d;                     % and small elements to the range D.m+1:D.n.
    if     D.i(l,d)==D.m, D.i(l,d)=D.i(D.m,d); D.i(D.m,d)=l;  % First, fix pointers...
    elseif D.i(D.m,d)==l, D.i(D.m,d)=D.i(l,d); D.i(l,d)=D.m;
    else   [D.i(l,d),D.i(D.m,d)]=Swap(D.i(l,d),D.i(D.m,d));  end
    if     D.k(l,d)==D.m, D.k(l,d)=D.k(D.m,d); D.k(D.m,d)=l;
    elseif D.k(D.m,d)==l, D.k(D.m,d)=D.k(l,d); D.k(l,d)=D.m;
    else   [D.k(l,d),D.k(D.m,d)]=Swap(D.k(l,d),D.k(D.m,d));  end
    D.i(D.k(l,d),d)=l; D.i(D.k(D.m,d),d)=D.m;     
	D.k(D.i(l,d),d)=l; D.k(D.i(D.m,d),d)=D.m;
  end                                             
  [D.P(l),  D.P(D.m)  ]=Swap(D.P(l),  D.P(D.m));    % ... then swap elements l and D.m.
  [D.j(l,:),D.j(D.m,:)]=Swap(D.j(l,:),D.j(D.m,:)); 
  [D.v(l,:),D.v(D.m,:)]=Swap(D.v(l,:),D.v(D.m,:)); 
  [D.f(l,:),D.f(D.m,:)]=Swap(D.f(l,:),D.f(D.m,:));
  [D.u(l,:),D.u(D.m,:)]=Swap(D.u(l,:),D.u(D.m,:));
  [D.w(l,:),D.w(D.m,:)]=Swap(D.w(l,:),D.w(D.m,:));
  D.m=D.m-1;
end, l=l-1; end   
D.f(D.m+1:D.n,1)=zeros(D.n-D.m,1);  % Next, identify the neighbors to the big elements,
for l=2:D.m, for d=1:G.d            % and create entries for them if necessary.
  if D.i(l,d)==1, D=Create(D,G,D.j(l,:)-G.Y(d,:)); else, D.f(D.i(l,d),1)=1; end
  if D.k(l,d)==1, D=Create(D,G,D.j(l,:)+G.Y(d,:)); else, D.f(D.k(l,d),1)=1; end
  for e=d:G.d,                      % (the following computes the corner entries)
    if D.i(D.i(l,e),d)==1, D=Create(D,G,D.j(l,:)-G.Y(d,:)-G.Y(e,:));
	  else, D.f(D.i(D.i(l,e),d),1)=1; end
    if D.i(D.k(l,e),d)==1, D=Create(D,G,D.j(l,:)-G.Y(d,:)+G.Y(e,:));
	  else, D.f(D.i(D.k(l,e),d),1)=1; end
    if D.k(D.i(l,e),d)==1, D=Create(D,G,D.j(l,:)+G.Y(d,:)-G.Y(e,:));
	  else, D.f(D.k(D.i(l,e),d),1)=1; end
    if D.k(D.k(l,e),d)==1, D=Create(D,G,D.j(l,:)+G.Y(d,:)+G.Y(e,:));
	  else, D.f(D.k(D.k(l,e),d),1)=1; end
  end
end, end  

l=D.m+1; while l<=D.n, if D.f(l,1)~=1,   % Remove small elements which do not neighbor
  for d=1:G.d;                           % the big elements.  First, fix pointers...
    D.i(D.k(l,d),d)=1; if l<D.n, D.i(D.k(D.n,d),d)=l; end
    D.k(D.i(l,d),d)=1; if l<D.n, D.k(D.i(D.n,d),d)=l; end
    if D.j(D.n,:)-G.Y(d,:)==D.j(l,:), D.i(l,d)=1; else, D.i(l,d)=D.i(D.n,d); end
    if D.j(D.n,:)+G.Y(d,:)==D.j(l,:), D.k(l,d)=1; else, D.k(l,d)=D.k(D.n,d); end
  end                                            
  D.j(l,:)=D.j(D.n,:); D.P(l)=D.P(D.n); D.v(l,:)=D.v(D.n,:);     % then replace element
  D.u(l,:)=D.u(D.n,:); D.w(l,:)=D.w(D.n,:); D.f(l,1)=D.f(D.n,1); % l with element D.n.
  D.n=D.n-1;
else, l=l+1; end, end
for l=2:D.n, D.P(l)=max(D.P(l),0); end, D.P(1)=0; D.i(1,1:G.d)=1; D.k(1,1:G.d)=1;
D.P(1:D.n,1)=D.P(1:D.n,1)/sum(D.P(1:D.m,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,a]=Swap(a,b); end  % Swap two elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=Create(D,G,j);  % Create a new element in D
D.n=D.n+1; D.P(D.n)=0;  D.j(D.n,:)=j; [D]=Initialize_vuwik(D,G,D.n); D.f(D.n,1)=1;
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=G_RHS_P(D,P,G)
D.f(1:D.n,1:G.d)=0; P(1)=0; K(1:D.n,1)=0;
% -------------------------------- NONCONSERVATIVE FORM -------------------------------
% for l=2:D.n, for d=1:G.d, K(l,1)=K(l,1)-D.w(D.i(l,d),d)*(P(l)-P(D.i(l,d)))/G.dx ...
%							             -D.u(l,d)*(P(D.k(l,d))-P(l))/G.dx;    end, end
% --------------------------------- CONSERVATIVE FORM ---------------------------------                                         
for l=2:D.n, for d=1:G.d, D.f(l,d)=D.w(l,d)*P(l)+D.u(l,d)*P(D.k(l,d)); end, end
%------------------------- (keep one of the above 2 sections) -------------------------
for d=1:G.d, for l=2:D.n, i=D.i(l,d); if l<=D.m | (i>1 & i<=D.m),
  F=G.dt*(P(l)-P(i))/(2*G.dx);      
  for e=1:G.d, if e~=d,
    D.f(l,e)=D.f(l,e)-D.w(l,e)*D.w(i,d)*F; j=D.i(l,e);   % Compute corner transport
    D.f(j,e)=D.f(j,e)-D.u(j,e)*D.w(i,d)*F;               % upwind (CTU) flux terms.
    D.f(i,e)=D.f(i,e)-D.w(i,e)*D.u(i,d)*F; p=D.i(i,e);
    D.f(p,e)=D.f(p,e)-D.u(p,e)*D.u(i,d)*F;
   end, end
  if D.v(i,d)>0, th=(P(i)-P(D.i(i,d)))/(P(l)-P(i));            % Compute second-order
  else,          th=(P(D.k(l,d))-P(l))/(P(l)-P(i)); end,       % correction flux term.
  t=abs(D.v(i,d)); D.f(i,d)=D.f(i,d)+t*(G.dx/G.dt-t)*F*MC(th); % Flux: use MC or VL.
end, end, end
for l=2:D.n, for d=1:G.d, K(l,1)=K(l,1)-(D.f(l,d)-D.f(D.i(l,d),d))/G.dx; end, end
end % function RHS_P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi]=MC(th), phi=max(0,min([(1+th)/2 2 2*th]));   end         % Flux limiters
function [phi]=VL(th), phi=min((th+abs(th))/(1+abs(th)),0); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     HGBEES_1 FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = H1_Initialize_D(G)
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
    D.m=numEntries(D.P); D.n=D.m; D = H1_Initialize_vuw(D,G,2); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=H1_Initialize_vuw(D,G,b)
    for l=b:D.n
        current_key = D.keys(l); current_state = double(key_conversion(current_key,G));
        x=G.dx.*current_state; xh=G.dx/2;
        v1=G.sigma*(x(2)-(x(1)+xh)); v2=-(x(2)+xh)-x(1)*x(3); v3=-G.b*(x(3)+xh)+x(1)*x(2)-G.b*G.r; 
        D.v(current_key) = {[v1 v2 v3]};
        D.u(current_key)={[min(v1, 0) min(v2, 0) min(v3, 0)]};
        D.w(current_key)={[max(v1, 0) max(v2, 0) max(v3, 0)]};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = state_conversion(state,G)
    state = int32(state); state = uint64(state + G.offset32);
    key=uint64(0);
    for i=1:G.N_data; key=key+state(i)*G.fac^(i-1); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = key_conversion(key,G)
    for i=G.N_data:-1:1
        state(i)=idivide(key,G.fac^(i-1),'floor');
        key=key-state(i)*G.fac^(i-1);
    end
    state=int64(state)-G.offset64;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = H1_Modify_pointset(D,G) 
    D.keys = keys(D.P); D.m = numEntries(D.P);
    for l=2:D.m    % Check/Create Neighbors of Big Cells
        if(D.P(D.keys(l))>=G.thresh)
            og_state = key_conversion(D.keys(l),G);
    
            %   Faces (6)
            new_key = state_conversion([og_state(1) + 1, og_state(2), og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) - 1, og_state(2), og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1), og_state(2) + 1, og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1), og_state(2) - 1, og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
            new_key = state_conversion([og_state(1), og_state(2), og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
            new_key = state_conversion([og_state(1), og_state(2), og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
    
            %   Edges (12)
            new_key = state_conversion([og_state(1) + 1, og_state(2) + 1, og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) + 1, og_state(2) - 1, og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
            new_key = state_conversion([og_state(1) - 1, og_state(2) + 1, og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end   
            new_key = state_conversion([og_state(1) - 1, og_state(2) - 1, og_state(3)],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) + 1, og_state(2), og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) + 1, og_state(2), og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) - 1, og_state(2), og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) - 1, og_state(2), og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1), og_state(2) + 1, og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1), og_state(2) + 1, og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1), og_state(2) - 1, og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1), og_state(2) - 1, og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            
            %   Corners (8)
            new_key = state_conversion([og_state(1) + 1, og_state(2) + 1, og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) + 1, og_state(2) + 1, og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) + 1, og_state(2) - 1, og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) - 1, og_state(2) + 1, og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) + 1, og_state(2) - 1, og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) - 1, og_state(2) - 1, og_state(3) + 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) - 1, og_state(2) + 1, og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = state_conversion([og_state(1) - 1, og_state(2) - 1, og_state(3) - 1],G); if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P);
    D = H1_Initialize_vuw(D,G,D.m+1); 
    for l=2:D.m                     % Remove Small Elements
        if(D.values(l) < G.thresh)&&(H1_no_neighbors(D,G,l))
            D.P(D.keys(l)) = []; D.v(D.keys(l)) = []; D.u(D.keys(l)) = []; D.w(D.keys(l)) = [];
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P);    
    D.P(D.keys) = max(D.P(D.keys), 0); D.values = values(D.P); 
    prob_sum = sum(values(D.P)); D.P(D.keys) = D.values./prob_sum; D.values = values(D.P);
end                                                                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=H1_RHS_P(D,G)
    D.n = numEntries(D.P); D.keys = keys(D.P); D.values = values(D.P); K(1:D.n,1)=0; 
    D.f(D.keys) = {zeros(1,G.d)};
    for l=2:D.n
        l_key = D.keys(l); l_state = key_conversion(l_key,G);
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
            l_key = D.keys(l); l_state = key_conversion(l_key,G);
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
        l_key = D.keys(l); l_state = key_conversion(l_key,G);
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
function bool=H1_no_neighbors(D,G,l)
    bool = 1;
    current_key = D.keys(l); og_state = key_conversion(current_key,G);
    % Faces (6)
    new_key = state_conversion([og_state(1)+1 og_state(2) og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2) og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)+1 og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)-1 og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = state_conversion([og_state(1) og_state(2) og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = state_conversion([og_state(1) og_state(2) og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    
    % Edges (12)
    new_key = state_conversion([og_state(1)+1 og_state(2)+1 og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)-1 og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = state_conversion([og_state(1)-1 og_state(2)+1 og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end   
    new_key = state_conversion([og_state(1)-1 og_state(2)-1 og_state(3)],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2) og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2) og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2) og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2) og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)+1 og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)+1 og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)-1 og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1) og_state(2)-1 og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    %   Corners (8)
    new_key = state_conversion([og_state(1)+1 og_state(2)+1 og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)+1 og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)-1 og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)+1 og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)+1 og_state(2)-1 og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)-1 og_state(3)+1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)+1 og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = state_conversion([og_state(1)-1 og_state(2)-1 og_state(3)-1],G); if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     HGBEES_2 FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = H2_Initialize_D(G)
    D.P = dictionary(); D.P({-10000*ones(1,G.d)}) = 0;
    D.v = dictionary(); D.v({-10000*ones(1,G.d)}) = {zeros(1,G.d)};
    D.u = D.v; D.w = D.v; D.f = D.v;
    D.neighbors = dictionary();

    for i=round((G.start(1)-2)/G.dx):round((G.start(1)+2)/G.dx)
        for j=round((G.start(2)-2)/G.dx):round((G.start(2)+2)/G.dx)
            for k=round((G.start(3)-2)/G.dx):round((G.start(3)+2)/G.dx)
                key = {[i j k]};
                x=(i*G.dx - G.start(1))^2+(j*G.dx - G.start(2))^2+(k*G.dx - G.start(3))^2;
                D.P(key) = exp(-4.*x/2.); 
            end
        end
    end
    D.keys = keys(D.P); D.m=numEntries(D.P); D.n=D.m; D = H2_Initialize_vuw(D,G,2); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D]=H2_Initialize_vuw(D,G,b)
    for l=b:D.n
        current_key = D.keys(l); current_state = current_key{1};
        x=G.dx.*(current_state); xh=G.dx/2;
        v1=G.sigma*(x(2)-(x(1)+xh)); v2=-(x(2)+xh)-x(1)*x(3); v3=-G.b*(x(3)+xh)+x(1)*x(2)-G.b*G.r; 
        D.v(current_key) = {[v1 v2 v3]};
        D.u(current_key)={[min(v1, 0) min(v2, 0) min(v3, 0)]};
        D.w(current_key)={[max(v1, 0) max(v2, 0) max(v3, 0)]};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D] = H2_Modify_pointset(D,G) 
    D.keys = keys(D.P); D.m = numEntries(D.P);
    for l=2:D.m    % Check/Create Neighbors of Big Cells
        if(D.P(D.keys(l))>=G.thresh)
            og_state = D.keys(l); og_state = og_state{1};
     
            %   Faces (6)
            new_key = {[og_state(1)+1 og_state(2) og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)-1 og_state(2) og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1) og_state(2)+1 og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1) og_state(2)-1 og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
            new_key = {[og_state(1) og_state(2) og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
            new_key = {[og_state(1) og_state(2) og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
            
            %   Edges (12)
            new_key = {[og_state(1)+1 og_state(2)+1 og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)+1 og_state(2)-1 og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end        
            new_key = {[og_state(1)-1 og_state(2)+1 og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end   
            new_key = {[og_state(1)-1 og_state(2)-1 og_state(3)]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)+1 og_state(2) og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)+1 og_state(2) og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)-1 og_state(2) og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)-1 og_state(2) og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1) og_state(2)+1 og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1) og_state(2)+1 og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1) og_state(2)-1 og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1) og_state(2)-1 og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end

            %   Corners (8)
            new_key = {[og_state(1)+1 og_state(2)+1 og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)+1 og_state(2)+1 og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)+1 og_state(2)-1 og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)-1 og_state(2)+1 og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)+1 og_state(2)-1 og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)-1 og_state(2)-1 og_state(3)+1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)-1 og_state(2)+1 og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
            new_key = {[og_state(1)-1 og_state(2)-1 og_state(3)-1]}; if(~isKey(D.P, new_key)), D.P(new_key) = 0; end
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P); 
    D = H2_Initialize_vuw(D,G,D.m+1);

    for l=2:D.m                     % Remove Small Elements which do not neighbor big elements
        if(D.P(D.keys(l)) < G.thresh)&&(H2_no_neighbors(D,G,l))
            D.P(D.keys(l)) = []; D.v(D.keys(l)) = []; D.u(D.keys(l)) = []; D.w(D.keys(l)) = [];
        end
    end
    D.n = numEntries(D.P); D.keys = keys(D.P);  
    D.P(D.keys) = max(D.P(D.keys), 0); prob_sum = sum(D.P(D.keys)); D.P(D.keys) = D.P(D.keys)./prob_sum;  
end                                                                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=H2_RHS_P(D,G)
    D.n = numEntries(D.P); D.keys = keys(D.P); K(1:D.n,1)=0;
    D.f(D.keys) = {zeros(1,G.d)};
    for l=2:D.n
        l_key = D.keys(l); l_state = l_key{1};
        f_l = D.f(l_key); f_l = f_l{1}; u_l = D.u(l_key); u_l = u_l{1}; w_l = D.w(l_key); w_l = w_l{1};  
        for d=1:G.d
            k_state = l_state; k_state(d) = k_state(d)+1; k_key = {k_state};
            if(~isKey(D.P,k_key)),k_key = {-10000*ones(1,G.d)};end
            f_l(d) = w_l(d)*D.P(l_key) + u_l(d)*D.P(k_key);
        end
        D.f(l_key) = {f_l};
    end

    for d=1:G.d
        for l=2:D.n
            l_key = D.keys(l); l_state = l_key{1};
            f_l = D.f(l_key); f_l = f_l{1}; w_l = D.w(l_key); w_l = w_l{1}; 
            i_state = l_state; i_state(d) = i_state(d)-1; i_key = {i_state};
            if(~isKey(D.P,i_key)),i_key = {-10000*ones(1,G.d)};end
            f_i = D.f(i_key); f_i = f_i{1}; w_i = D.w(i_key); w_i = w_i{1}; u_i = D.u(i_key); u_i = u_i{1};
            v_i = D.v(i_key); v_i = v_i{1};
            if (D.P(l_key)>=G.thresh)||(D.P(i_key)>=G.thresh)
                F=G.dt*(D.P(l_key)-D.P(i_key))/(2*G.dx);
                for e=1:G.d
                    if e~=d
                        j_state = l_state; j_state(e) = j_state(e)-1; j_key = {j_state};
                        if(~isKey(D.P,j_key)),j_key = {-10000*ones(1,G.d)};end
                        f_j = D.f(j_key); f_j = f_j{1}; u_j = D.u(j_key); u_j = u_j{1};
                        p_state = i_state; p_state(e) = p_state(e)-1; p_key = {p_state};
                        if(~isKey(D.P,p_key)),p_key = {-10000*ones(1,G.d)};end
                        f_p = D.f(p_key); f_p = f_p{1}; u_p = D.u(p_key); u_p = u_p{1};

                        f_l(e) = f_l(e)-w_l(e)*w_i(d)*F; D.f(l_key) = {f_l};
                        f_j(e) = f_j(e)-u_j(e)*w_i(d)*F; D.f(j_key) = {f_j};
                        f_i(e) = f_i(e)-w_i(e)*u_i(d)*F; D.f(i_key) = {f_i};
                        f_p(e) = f_p(e)-u_p(e)*u_i(d)*F; D.f(p_key) = {f_p};
                    end
                end                
               
                i_i_state = i_state; i_i_state(d) = i_i_state(d)-1; i_i_key = {i_i_state};
                if(~isKey(D.P,i_i_key)),i_i_key = {-10000*ones(1,G.d)};end
                k_state = l_state; k_state(d) = k_state(d)+1; k_key = {k_state};
                if(~isKey(D.P,k_key)),k_key = {-10000*ones(1,G.d)};end
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
        l_key = D.keys(l); l_state = l_key{1};
        f_l = D.f(l_key); f_l = f_l{1};
        for d=1:G.d
            i_state = l_state; i_state(d) = i_state(d)-1; i_key = {i_state};
            if(~isKey(D.P,i_key)),i_key = {-10000*ones(1,G.d)};end
            f_i = D.f(i_key); f_i = f_i{1};
            K(l,1)=K(l,1)-(f_l(d)-f_i(d))/G.dx;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool=H2_no_neighbors(D,G,l)
    bool = 1;
    current_key = D.keys(l); og_state = current_key{1};
    % Faces (6)
    new_key = {[og_state(1)+1 og_state(2) og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)-1 og_state(2) og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1) og_state(2)+1 og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1) og_state(2)-1 og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = {[og_state(1) og_state(2) og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = {[og_state(1) og_state(2) og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    
    % Edges (12)
    new_key = {[og_state(1)+1 og_state(2)+1 og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)+1 og_state(2)-1 og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end        
    new_key = {[og_state(1)-1 og_state(2)+1 og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end   
    new_key = {[og_state(1)-1 og_state(2)-1 og_state(3)]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)+1 og_state(2) og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)+1 og_state(2) og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)-1 og_state(2) og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)-1 og_state(2) og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1) og_state(2)+1 og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1) og_state(2)+1 og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1) og_state(2)-1 og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1) og_state(2)-1 og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end

    %   Corners (8)
    new_key = {[og_state(1)+1 og_state(2)+1 og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)+1 og_state(2)+1 og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)+1 og_state(2)-1 og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)-1 og_state(2)+1 og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)+1 og_state(2)-1 og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)-1 og_state(2)-1 og_state(3)+1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)-1 og_state(2)+1 og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
    new_key = {[og_state(1)-1 og_state(2)-1 og_state(3)-1]}; if(isKey(D.P, new_key)), if(D.P(new_key)>=G.thresh), bool = 0; return; end, end
end