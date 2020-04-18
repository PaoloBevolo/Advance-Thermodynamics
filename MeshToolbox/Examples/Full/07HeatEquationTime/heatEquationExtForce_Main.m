clear
TRef=1; %"reference" temperature
%% Geometry definition & boundary conditions
Reg=regions.rect('mu',237,'rho',2700*896.9)-regions.circle([0,0],[0.1,0.1],64);
figure;
Reg.draw('e');

Reg.Borders(1).Bc([1,3])=boundaries.neumann(0);
Reg.Borders(2).Bc(:)=boundaries.dirichlet(20);
figure;
Reg.draw('bc');
%% Mesh
Me=mesh2D(Reg,0.001);
figure;
Me.draw('d')
%% Stifness matrix 
[D, bBoundary, bExt]=heatEquationExtForce_BuildStiff(Me,@(x,y)ones(size(x)));

%% Stationary solution: 
InitialForce=0;
uStationary0=D\(bBoundary+bExt*InitialForce);
%solution plot
uu=Me.copyToAllNodes(uStationary0);
figure; 
Me.draw(uu);
ylabel(colorbar(),'Temperature [^\circC]');
%% time evolution: definition of external force
%edge 4 temperature
Tend=4000;
f=@(t)InitialForce+(t<1000).*(30*t)+(t>=1000).*30000;
figure;fplot(f,[0,3000]);
xlabel('Time [s]');
ylabel('fExt');
grid on;

M=buildMass(Me);
Dofs=Me.Nodes.Dof>0;
%% implicit Euler method
disp('Implicit Euler method');
u=uStationary0;
uu=Me.copyToAllNodes(uStationary0);
dt=100;
figure;
S=(M+D*dt);
tic
for k=1:Tend/dt,
    t=k*dt;
    %original: u=(M+D*dt)\(M*u+b*dt*T4);
    u=S\(M*u+dt*(bBoundary+bExt*f(t)));
    uu(Dofs)=u;
    hold off;    
    Me.draw(uu,'hidemesh');
    title(['t = ' num2str(t) 's']);   
    zlim([0 20])
    caxis([0 20])
    xlabel('x dir [m]');
    ylabel('y dir [m]');
    ylabel(colorbar(),'Temperature [^\circC]');
    drawnow();
end
toc

%% Crank-Nicholson method (
dt=100;
u=uStationary0; 
uu=Me.copyToAllNodes(uStationary0);
B=D*dt/2;  
CN1=M+B;  
CN2=M-B; 
Fold=f(0);
for k=1:Tend/dt,	
    t=k*dt;	
    F=f(t);   
    u=CN1\(CN2*u+dt*(bBoundary+bExt*(F+Fold)/2));	
    uu(Dofs)=u;	
    Me.draw(uu,'hidemesh');
    hold off; 
    zlim([0 20]);	
    caxis([0,20]);
    title(['t=',num2str(t) 's']);
    drawnow();   	
    Fold=F; 
end

%% Explicit Euler method
dtMax=2/eigs(D,M,1,'lm');
dt=dtMax/1.5;%more safe
u=uStationary0;
uu=Me.copyToAllNodes(uStationary0);
Ddt=D*dt;		
for k=1:Tend/dt,    
    t=k*dt;    
    F=f(t-dt); %previous step! 	 
    u=u-M\(Ddt*u-dt*(bBoundary+bExt*F));	 
    if rem (k,100)==0 %plot once every 100 step
        uu(Dofs)=u;	 
        Me.draw(uu,'hidemesh');	 
        hold off; 	 
        zlim([0 20]);	 
        caxis([0,20]);	 
        title(['t=',num2str(t) 's']); 
        drawnow(); 
    end
end

%% Ode45
u0 = uStationary0;
fode = @(t,u)-D*u+(bBoundary+bExt*f(t));
odepar = odeset('Mass',M); %Mass matrix
[t,U] = ode45(fode, [0,Tend],u0,odepar);
size(U) % 32069 x 501 ->501 time steps
figure;
for k = 1:100:length(t)	
    uu(Dofs) = U(k,:);	
    Me.draw(uu,'hidemesh');	
    hold off;
    zlim([0 20]);	
    caxis([0 20]);	%  view([0 90]);	
    title(['t= ' num2str(t(k)) 's']);	
    drawnow();  
end





