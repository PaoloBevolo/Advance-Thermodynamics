%Dirichlet non homogeneous B.C.s
Reg=regions.rect('mu',1)+[0.5 0.5];
Reg.draw('edge')
f=@(x,y)-4* ones(size(x));
Reg.Borders.Bc(3)=boundaries.dirichlet(1);
Reg.Borders.Bc(4)=boundaries.dirichlet(@(x,y)x);
Reg.Borders.Bc(2)=boundaries.dirichlet(@(x,y)x);
Me=mesh2D(Reg, 0.05);
[A,b]=dirichletNonHomo_BuildStiff(Me,f);
u=pcg(A,b,1e-3,100);

figure;
Me.draw()
%uu=zeros(size(Me.Nodes.Dof));
%uu(Me.Nodes.Dof<0)=Me.BC.DirichletNodes(:,2);
%uu(Me.Nodes.Dof>0)=u; %impongo la soluzione nei nodi interni
uu=Me.copyToAllNodes(u);
figure;
subplot(2,1,1);
Me.draw(uu,'offset');
subplot(2,1,2);
Me.draw(uu,'contour');