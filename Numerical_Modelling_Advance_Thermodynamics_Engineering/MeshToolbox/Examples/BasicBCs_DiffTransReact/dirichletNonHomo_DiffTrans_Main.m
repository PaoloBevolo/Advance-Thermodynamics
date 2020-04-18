%Example of solution with Homogeneous Dirichlet B.C.s.
%Transport-diffusion-reaction problem

Reg=regions.rect('mu',1,'beta',25*[1,-1]);
Reg.Borders.Bc(:)=boundaries.dirichlet(2);
f=@(x,y)zeros(size(x))+4;   %external force
Me=mesh2D(Reg,0.0001);

%Stiffness matrix construction
[DT,b]=dirichletNonHomo_DiffTrans_BuildStiff(Me,f); 
uu=Me.copyToAllNodes(DT\b);
figure;
Me.draw(uu,'hidemesh');
shading ('interp');
xlabel ('x axis')
ylabel ('y axis')
view ([135 25]);

