Reg(1)=regions.rect([0 0.5],[10,1],'mu',1)-regions.rect([0 0.5],[3,0.3]);
Reg(2)=regions.rect([0 -1],[8,2],'mu',2)-regions.circle([0 -1],0.5);
Reg(3)=regions.rect([0 0.5],[3,0.3],'mu',1e-3);
figure;
Reg.draw();
Reg(1).Borders(1).insertNode(1,[4 0;-4 0]);
figure;
Reg.draw('nodes');
Reg(1).Borders(1).Bc([1 3 5])=boundaries.neumann(0);
Reg(1).Borders(1).Bc(2)=boundaries.continuity();
Reg(1).Borders(1).Bc([6 4])=boundaries.periodic(@(x,y)[-x,y]);
Reg(1).Borders(2).Bc(:)=boundaries.continuity();
Reg(2).Borders(1).Bc([1 2 4])=boundaries.neumann(0);
Reg(2).Borders(1).Bc(3)=boundaries.continuity();
Reg(3).Borders(1).Bc(:)=boundaries.continuity();
figure;
Reg.draw('bc');
xyc=[0 0.5]; %centro del rettangolo
wh=[4,0.6]; %larghezza ed altezza
Me=mesh2D(Reg,0.01);
figure;
Me.draw('d');
%%
f=@(x,y)1*InRect([x,y],[-3,0.5],[1,0.5])-1*InRect([x,y],[-2.5,-0.5],[1,0.5]);
[A,b]=periodic_BuildStiff(Me,f);
%uu=zeros(size(Me.Coordinates,1),1);
%uu(Me.UnknownNodes>0)=  pcg(A,b,1e-3,1000);
uu=Me.copyToAllNodes(A\b);
figure;
Me.draw(uu,'offset');
title('Periodic conditions');
