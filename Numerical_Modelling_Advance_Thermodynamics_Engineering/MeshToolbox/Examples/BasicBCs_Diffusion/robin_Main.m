%Example of solution with Robin B.C.
Reg=regions.rect('mu',1);
%tending to g/h with exponential [exp(-hx)] behavior
Reg.Borders.Bc(1)=boundaries.robin(@(x,y)2,@(x,y)-2);
f=@(x,y)zeros(size(x));
Me=mesh2D(Reg,0.00001);
[D,b]=robin_BuildStiff(Me,f);
u=D\b;
uu=Me.copyToAllNodes(u);
figure;
Me.draw(uu,'hidemesh');
shading interp;
