%Darcy problem solved in stationary conditions
grid2D.setSteps(0.0001);
Reg=regions.rect([-1, 3.5],[2, 1],'mu',1)+...
    regions.sector([0, 0],4,3,[0, 90],16)+...
    regions.rect([3.5, -1],[1, 2]);
figure;
Reg.draw('e')
axis('equal');

Reg.Borders.Bc(:)=boundaries.neumann(0);
Reg.Borders.Bc(36)=boundaries.neumann(4);
Reg.Borders.Bc(17)=boundaries.neumann(-4);
figure;
Reg.draw('bc');

Me=mesh2D(Reg,0.02);
DirNode=Me.findClosestNode([-2,4]);
Me=Me.forceDirichlet(DirNode, 0);
figure;
Me.draw('d');
[A,b]=neumannNonHomo_BuildStiff(Me,@(x,y)zeros(size(x)));

uu=Me.copyToAllNodes(A\b);
figure
Me.draw(uu)
title('Pressure profile');

[ux,uy]=Me.gradient(uu);
figure; 
subplot(2,1,1);
Me.draw(ux,'hidemesh');
title('x-direction derivative');
subplot(2,1,2);
Me.draw(uy,'hidemesh');
title('y-direction derivative');

figure;
quiver(Me.Nodes.X,Me.Nodes.Y,-ux,-uy);
title('Velocity');
xlabel('x-dir');
ylabel('y-dir');

