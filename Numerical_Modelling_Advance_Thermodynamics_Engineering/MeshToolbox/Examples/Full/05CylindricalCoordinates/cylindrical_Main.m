%Example of solution of the homogeneous Dirichlet problem with cylindrical
%coordinates
grid2D.setSteps(1e-5);
Reg=regions.rectN([-1,0],[1,0.5],'mu',5,'beta',[1,0])-regions.circle([0,0],0.1,16);
figure;
subplot(2,1,1);
Reg.draw('e');
axis equal
Reg.Borders.Bc([1,10])=boundaries.neumann(0);
Reg.Borders.Bc(2:9)=boundaries.dirichlet(20);
Reg.Borders.Bc(11:13)=boundaries.dirichlet(10);
subplot(2,1,2);
Reg.draw('bc');
Me=mesh2D(Reg,0.002);
%Cylindrical coordinates
[D,b]=cylindrical_BuildStiff(Me); %Assemblo la matrice di rigidezza
uuCylinder= Me.copyToAllNodes(D\b);

%Rectangular coordinates
[D,b]=dirichletNonHomo_DiffTrans_BuildStiff(Me,@(x,y)zeros(size(x))); %Assemblo la matrice di rigidezza
uuRectangular=Me.copyToAllNodes(D\b);
%
figure;
subplot(2,1,1);
Me.draw(uuCylinder,'c');
axis equal
xlim([-1, 1]);ylim([0, 0.5]);caxis([10 20]);
title('Cylindrical coordinates');

subplot(2,1,2);
Me.draw( uuRectangular,'c');
axis equal
xlim([-1, 1]);ylim([0, 0.5]);caxis([10 20]);
title('Cartesian coordinates');
