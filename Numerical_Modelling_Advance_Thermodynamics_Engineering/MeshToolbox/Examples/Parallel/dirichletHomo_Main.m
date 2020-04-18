%Solution with Homogeneous Dirichlet conditions 
Sh=shapes.rect('c',1);
Me=mesh2D(Sh,0.005);
%%
f=@(x,y)-4*ones(size(x));                           %External force
tic
[D,b]=dirichletHomo_BuildStiff(Me,f);   %Stiffness matrix
toc
%%
figure; spy(D)
title('Stiffness matrix pattern');
u=pcg(D,b,1e-6, 1000);                        %linear system solution
% uu=zeros(size(Me.Coordinates,1),1);   %alloco il vettore colonna soluzione
% uu(Me.UnknownNodes>0)=u;              %impongo la soluzione nei nodi interni
%%
uu=Me.copyToAllNodes(u);
figure;
subplot(2,1,1);
Me.draw(uu,'hidemesh');
subplot(2,1,2);
Me.draw(uu,'contour');
