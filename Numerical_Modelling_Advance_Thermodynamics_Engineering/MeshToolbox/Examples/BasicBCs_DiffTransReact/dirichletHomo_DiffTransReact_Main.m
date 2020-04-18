%Diffusion + reaction + transport problem with homogeneous Dirichlet
%conditions
Reg=regions.rect('mu',1,'sigma',5,'beta',25*[1,-1]);
f=@(x,y)4*ones(size(x));   %definisco la forzante
Me=mesh2D(Reg,0.02);
[D,T,R,b]=dirichletHomo_DiffTransReact_BuildStiff(Me,f); %Assemblo la matrice di rigidezza
figure; 
%% Diffusion only
subplot(2,2,1);
u=D\b;                      %Linear system solution
uu=Me.copyToAllNodes(u);
Me.draw(uu);
shading interp
title('Diffusion only');
view ([135 25]);

%% Reaction only
subplot(2,2,3);
u=R\b;                      %Linear system solution
uu=Me.copyToAllNodes(u);
Me.draw(uu);
shading interp
view ([135 25]);
title('Reaction only');

%% Transport only
subplot(2,2,2);
u=T\b;                      %Linear system solution
uu=Me.copyToAllNodes(u);
Me.draw(uu);
shading interp
view ([135 25]);
title('Transport only');

%% Diffusion + reaction + transport
subplot(2,2,4);
u=(D+T+R)\b;                      %Linear system solution
uu=Me.copyToAllNodes(u);
Me.draw(uu);
shading interp
view ([135 25]);
title('All together');