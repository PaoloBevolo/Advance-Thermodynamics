%% Solution with Homogeneous Dirichlet conditions 
%Create a rectangle centered in the origin, with unitary side, 
%and the parameter 'mu' set to 1
Reg=regions.rect('mu',1);
%open an new figure
figure;
%divide the figure in 2 rows and 2 columns; place the axes in position 1
subplot(2,2,1);
%Draw the geometry
Reg.draw(); 
title('Sh.draw()');
subplot(2,2,2);
Reg.draw('n');
title('Sh.draw(''n'')');
subplot(2,2,3);
Reg.draw('e');
title('Sh.draw(''e'')');
subplot(2,2,4);
Reg.draw('bc');
title('Sh.draw(''bc'')');
%%
Me=mesh2D(Reg,0.005);
figure;
subplot(2,3,1);
Me.draw();
title('Me.draw();');
subplot(2,3,2);
Me.draw('n');
title('Me.draw(''n'');');
subplot(2,3,3);
Me.draw('d');
title('Me.draw(''d'');');
subplot(2,3,4);
Me.draw('i');
title('Me.draw(''i'');');
subplot(2,3,5);
Me.draw('t');
title('Me.draw(''t'');');
subplot(2,3,6);
Me.draw('q');
title('Me.draw(''q'');');
%%
f=@(x,y)-4*ones(size(x));               %External force
[D,b]=dirichletHomo_BuildStiff(Me,f);   %Stiffness matrix
figure; spy(D)
title('Stiffness matrix pattern');
u=pcg(D,b,1e-3);                        %linear system solution

uu=Me.copyToAllNodes(u);
figure;
subplot(2,1,1);
Me.draw(uu,'hidemesh');
subplot(2,1,2);
Me.draw(uu,'contour');

%% non uniform elasticity
Reg(1)=regions.rect();
Reg(2)=regions.rect()+[1,0];
Reg=Reg(1)+Reg(2);
Reg.addProperty('mu',@(x,y)(1+(x>0.5)));
Me=mesh2D(Reg,0.005);
[A,b]=dirichletHomo_BuildStiff(Me, @(x,y)-4*ones(size(x)));
u=pcg(A,b,1e-3,200);
% uu=zeros(size(Me.Coordinates,1),1);
% uu(Me.UnknownNodes<0)=0;
% uu(Me.UnknownNodes>0)=u;
uu=Me.copyToAllNodes(u);
figure;
Me.draw(uu,'offset');

%% convergence
N=100;
x=linspace(-0.5,1.5,N);
xy=[x.',zeros(N,1)];
figure; 
hold all;
Meshes=[0.02, 0.01, 0.005, 0.0025,0.00125];
zold=[];
Differences=zeros(numel(Meshes)-1,1);
for kIndex=1:numel(Meshes)
    maxArea=Meshes(kIndex);
    Me=mesh2D(Reg,maxArea);
    [A,b]=dirichletHomo_BuildStiff(Me, @(x,y)-4*ones(size(x)));
    u=pcg(A,b,1e-4,2000);
    uu=Me.copyToAllNodes(u);
    z=Me.interpolate(uu,xy);
    plot(x,z,'DisplayName',num2str(maxArea));
    if kIndex>1
        Differences(kIndex-1)=norm(zold-z);
    end
    zold=z;
    grid
end
legend('show');
xlabel('Position [m]')
ylabel('Displacement [m]')
disp(Differences);

%% two regions
%shape 1: square centered in the origin, 'mu' is set to 1
Reg(1)=regions.rect('mu',1);
%shape 2: square centered in the (1,0), 'c' is set to 1.5
Reg(2)=regions.rect('mu',1.5)+[1,0];
%Generate the mesh
Me=mesh2D(Reg,0.005);
%Build the stiffness matrix and the rhs. The force is -4 on the whole
%domain.
[A,b]=dirichletHomo_BuildStiff(Me, @(x,y)-4*ones(size(x)));
%Solve the linear system with pcg, 1e-3 relative tolerance, maximum 200
%iterations
u=pcg(A,b,1e-3,200);
%assign the solution vector to the nodes
uu=Me.copyToAllNodes(u);
%open a new figure
figure;
%draw the solution
Me.draw(uu);

