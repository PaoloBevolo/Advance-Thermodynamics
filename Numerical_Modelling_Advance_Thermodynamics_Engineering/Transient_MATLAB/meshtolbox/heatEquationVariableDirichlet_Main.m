clear
TRef=1; %"reference" temperature
%% Geometry definition & boundary conditions
Reg=regions.rect('mu',237,'rho',2700*896.9)-regions.circle([0,0],[0.2,0.2],64);
figure;
Reg.draw('e');

Reg.Borders(1).Bc([1,3])=boundaries.neumann(0);
Reg.Borders(1).Bc([2,4])=boundaries.dirichlet(20);
Reg.Borders(2).Bc(:)=boundaries.dirichlet(1);%assume 1 degree
figure;
Reg.draw('bc');
%% Mesh
Me=mesh2D(Reg,0.001);
%% Nodes on the circle, with Dirichlet B.C.
DirichletNodesCircle=Me.find(@(x,y)abs(x)<=0.2,'d');
DirichletNodesExtBorder=Me.find(@(x,y)abs(x)>0.2,'d');
%% Degrees of freedom
Dof=Me.Nodes.Dof>0;
%% Stifness matrix 
[D, bconst, bvar]=heatEquationVariableDirichlet_BuildStiff(Me);

%% Stationary solution: 
%bvar was calculated with Tref=1, now I have T0
T0=20; %initial temperature on the inner border
uStationary0=D\(bconst+bvar*T0);
%solution plot
uu=zeros(size(Me.Nodes.X));
uu(Dof>0)=uStationary0;
uu(DirichletNodesCircle)=T0;
uu(DirichletNodesExtBorder)=20;
figure; 
Me.draw( uu);
ylabel(colorbar(),'Temperature [^\circC]');
%% time evolution: definition of the function T
%edge 4 temperature
Tend=1000;
fDirichlet=@(t)20+(t<10)*8*t+80*(t>=10);
[M, mvar] = buildMassVariableDirichlet(Me);

%% implicit Euler method
disp('Implicit Euler method');
u=uStationary0;
dt=1;
figure;
A=(M+D*dt);
tic
TCircleOld=fDirichlet(0);
for k=1:Tend/dt,
    t=k*dt;
    TCircle=fDirichlet(t);
    %original: u=(M+D*dt)\(M*u+b*dt*T4);
    u=A\(M*u+dt*(bconst+bvar*TCircle)-mvar*(TCircle-TCircleOld));
    uu(Dof)=u;
    uu(DirichletNodesCircle)=TCircle;
    TCircleOld=TCircle;
    hold off;    
    if rem (k,10)==0
        Me.draw(uu,'hidemesh');
        zlim([0 100]);
        caxis([0 100]);
        title(['t = ' num2str(t) 's']);    
        xlabel('x dir [m]');
        ylabel('y dir [m]');
        ylabel(colorbar(),'Temperature [^\circC]');
        drawnow();
    end
end
toc
%% implicit Euler method, with pre-conditioning
disp('Implicit Euler method, with pre-conditioning');
u=uStationary0;
dt=1;
figure;
A=(M+D*dt);
R=ichol(A,struct('type','ict','droptol',1e-3,'shape','upper'));
tic
TCircleOld=fDirichlet(0);
for k=1:Tend/dt,
    t=k*dt;
    TCircle=fDirichlet(t);
    %original: u=(M+D*dt)\(M*u+b*dt*T4);
    %u=A\(M*u+dt*(bconst+bvar*TCircle)-mvar*(TCircle-TCircleOld));
    u=pcg(A,M*u+dt*(bconst+bvar*TCircle)-mvar*(TCircle-TCircleOld),1e-4,1000,R,R',u);
    uu(Dof)=u;
    uu(DirichletNodesCircle)=TCircle;
    TCircleOld=TCircle;
    hold off;    
    if rem (k,10)==0
        Me.draw(uu,'hidemesh');
        zlim([0 100]);
        caxis([0 100]);
        title(['t = ' num2str(t) 's']);    
        xlabel('x dir [m]');
        ylabel('y dir [m]');
        ylabel(colorbar(),'Temperature [^\circC]');
        drawnow();
    end
end
toc

%% ODE45
disp('ODE45; please wait: the results will be shown at the end of the ODE integration');
u0=zeros(size(b));
figure;
fode=@(t,u) (-D*u+b*fDirichlet(t));
odepar= odeset('Mass',M,'reltol',1e-3);
[t,U]=ode45(fode, [0,Tend],u0,odepar);
size(U) % 32069 x 501

for k=1:length(t)
    uu(Dof)=U(k,:);
    uu(NodesEdge4)=fDirichlet(t(k));
    hold off;    
    Me.draw(uu,'hidemesh');
    zlim([0 TMax]);
    caxis([0 TMax]);
    %  view([0 90]);
    title(['t= ' num2str(t(k)) 's']);
    xlabel('x dir [m]');
    ylabel('y dir [m]');
    ylabel(colorbar(),'Temperature [^\circC]');
    drawnow();  
end
figure;
plot(diff(t));
xlabel('Step');
ylabel('Integration step [s]');
grid on;

