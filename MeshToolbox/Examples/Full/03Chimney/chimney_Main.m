%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%STATIONARY CASE%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reg=[...
    regions.rect('mu',0.8,'name', 'bricks')-regions.circle([0,0],0.25);...
    regions.circle([0,0],0.25,'mu',260,'name','metal')-regions.circle([0,0],0.10)...
    ];
figure;
Reg.draw('edge');
%%

Reg(1).Borders(1).Bc([1,3])=boundaries.neumann(0);
Reg(1).Borders(1).Bc(2)=boundaries.dirichlet(10);
Reg(1).Borders(1).Bc(4)=boundaries.robin(10, 10*10);
Reg(1).Borders(2).Bc(:)=boundaries.continuity();
Reg(2).Borders(2).Bc(:)=boundaries.dirichlet(40);
Reg(2).Borders(1).Bc(:)=boundaries.continuity();
figure;
Reg.draw('bc');

%%
Me=mesh2D(Reg,0.03);
figure; 
Me.draw('dirichlet');
axis equal
%%
[D, bconst, bvar]=chimney_BuildStiff(Me);
b=bconst+bvar;
TT=Me.copyToAllNodes(D\b);
figure; 
Me.draw(TT);

%%
x=-0.5:0.01:0.5;
figure; 
hold on;
col='rgbkmc';
deltay=[0; 0.25; 0.5];
for k=1:length(deltay)
    y=zeros(1,length(x))+deltay(k);
    TDiag=Me.interpolate(TT, [x;y]');
    plot(x,TDiag, col(k))
end
legend(num2str(deltay))
xlabel('Diagonal [m]');
ylabel('Temperature [^\circC]');
grid

%%

[dTdx, dTdy]=Me.gradient(TT);
figure; 
subplot(2,1,1);
Me.draw(-dTdx,'hidemesh')
view([45 45]);
subplot(2,1,2);
Me.draw(-dTdy,'hidemesh')
view([45 45]);

figure;
quiver(Me.Nodes.X,Me.Nodes.Y,dTdx, dTdy);
xlabel('X [m]');
ylabel('Y [m]');

return;
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%CASO DINAMICO%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sh=[...
    regions.rect('c',0.8,'rho',1,'name', 'bricks')-regions.circle([0,0],0.25); ...
    regions.circle([0,0],0.25,'c',260,'rho',1,'name','metal')-regions.circle([0,0],0.10)...
    ];
figure; Sh.draw('edge');
%%

Sh(1).Borders(1).Bc([1,3])=boundaries.neumann(0);
Sh(1).Borders(1).Bc(2)=boundaries.dirichlet(10);
TInf=10;
Sh(1).Borders(1).Bc(4)=boundaries.robin(10, 10*TInf);
Sh(1).Borders(2).Bc(:)=boundaries.continuity();
Sh(2).Borders(2).Bc(:)=boundaries.dirichlet(1);
Sh(1).Borders(2).Bc(:)=boundaries.continuity();
figure; 
Sh.draw('bc');
axis equal
%%
Me=mesh2D(Sh, 0.03);
figure; 
Me.draw('dirichlet');
axis equal
NodesPipe=Me.find(@(x,y)x~=-0.5,'d');
NodesEdge=Me.find(@(x,y)x==-0.5,'a');

%%
[D,bconst, bvar]=chimney_BuildStiff(Me);
T=D\(bconst+bvar*10);
Dof=Me.Nodes.Dof>0;
TT=zeros(size(Dof));
TT(Dof)=T;
TT(NodesPipe)=10;
TT(NodesEdge)=10;

figure; 
Me.draw(TT,'hidemesh');

f=@(t)10+30*t/600*(t<600)+30*(t>=600);
M=buildMass(Me);
%% implicit Euler method
dt=10;
Tend=900;
TMax=40;
figure;
A=(M+D*dt);
for k=1:Tend/dt,
    t=k*dt;
    PipeTemperature=f(t);
    T=A\(M*T+dt*(bconst+bvar*PipeTemperature));
    TT(Dof)=T;
    TT(NodesPipe)=PipeTemperature;
    hold off;    
    Me.draw(TT,'hidemesh');
    zlim([0 TMax]);
    caxis([0 TMax]);
    %view ([0 90]);
    title(['t= ' num2str(t) 's']);    
    drawnow();
end