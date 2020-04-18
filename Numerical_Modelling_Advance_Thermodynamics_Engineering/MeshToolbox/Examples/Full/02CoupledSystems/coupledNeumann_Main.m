
Reg1=regions.rect('mu',1,'sigma',10);
Me1=mesh2D(Reg1,0.01);
figure;
Me1.draw();
[D1,R1,b1]=coupledNeumann_BuildStiff(Me1);
Reg2=Reg1;
Reg2.Borders.Bc([1,2,3])=boundaries.dirichlet(2);
Reg2.Borders.Bc(4)=boundaries.neumann(0);
Reg2=Reg2.addProperty('mu',2);
Me2=mesh2D(Reg2,0.01);
%Me2=mesh2D(Me1,Sh2);
[D2,R2,b2]=coupledNeumann_BuildStiff(Me2);


%% metodo diretto
Atotale=[D1+R1 -R1; -R2 D2+R2];
btotale=[b1;b2];
uu=Atotale\btotale;
uu1=uu(1:end/2);
uu2=uu(end/2+1:end);
figure;
Me1.draw(uu1);
hold on;
Me2.draw(uu2);

%% metodo iterativo

uu2=2*ones(size(Me1.Nodes.Dof));

toll=1e-3;
itermax=100;

iter=0;
err=1;
poscenter=Me2.findClosestNode([0,0]);
poscenterlocal=Me2.Nodes.Dof(poscenter);
uu1old=0;
while err>toll && iter<itermax
    uu1=(D1+R1)\(b1+R1*uu2);
    uu2=(D2+R2)\(b2+R2*uu1);
    iter=iter+1;
    err= abs(uu1(poscenterlocal)-uu1old);
    fprintf('Iteration number:\t%d\terror:\t%8.5f%%\n',iter, err*100);
    uu1old=uu1(poscenterlocal);
end
figure;
Me1.draw(uu1);
hold on;
Me2.draw(uu2);


