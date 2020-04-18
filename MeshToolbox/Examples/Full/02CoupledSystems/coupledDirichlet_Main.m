Reg1=regions.rect('mu',1,'sigma',10);
Me1=mesh2D(Reg1,0.01);
figure;
Me1.draw();
[D1,R1,b1]=coupledDirichlet_BuildStiff(Me1);
Reg2=Reg1;
Reg2.Borders.Bc(:)=boundaries.dirichlet(2);
Reg2=Reg2.addProperty('sigma',3);
Me2=mesh2D(Reg2,0.01);
[D2,R2,b2]=coupledDirichlet_BuildStiff(Me2);

%% metodo diretto
Atotale=[D1+R1 -R1; -R2 D2+R2];
btotale=[b1;b2];
u=Atotale\btotale;

Dof1=Me1.Nodes.Dof>0;
uu1(Dof1)=u(1:end/2);
uu1(~Dof1)=0;
Dof2=Me2.Nodes.Dof>0;
uu2(Dof2)=u(end/2+1:end);
uu2(~Dof2)=2;
figure;
Me1.draw(uu1);
hold on;
Me2.draw(uu2);
title('Direct solution');
%% metodo iterativo
u2=2*ones(sum(Dof1),1);

toll=1e-3;
itermax=100;

iter=0;
err=1;
poscenter=Me2.findClosestNode([0,0]);
poscenterlocal=Me2.Nodes.Dof(poscenter);
u1old=0;
while err>toll && iter<=itermax
    
    u1=(D1+R1)\(b1+R1*u2);
    u2=(D2+R2)\(b2+R2*u1);
    iter=iter+1;
    err= abs(u1(poscenterlocal)-u1old);
    fprintf('Iteration number:\t%d\terror:\t%8.5f%%\n',iter, err*100);
    u1old=u1(poscenterlocal);
end
if err<=toll
    disp('Done!');
else
    disp('Requested tolerance not reached');
end

uu1(Dof1)=u1;
uu1(~Dof1)=0;
uu2(Dof2)=u2;
uu2(~Dof2)=2;

figure;
Me1.draw(uu1);
hold on;
Me2.draw(uu2);
title('Iterative solution');

