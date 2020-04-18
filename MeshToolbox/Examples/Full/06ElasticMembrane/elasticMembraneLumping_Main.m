%script per la risoluzione del problema temporale della membrana elastica

Reg = regions.rect('mu', 2, 'rho', 1);
Me = mesh2D(Reg,0.001);
f = @(x,y) zeros(size(x))-4;           %definisce la densita' di forza
[A, b] = dirichletHomo_BuildStiff(Me, f);   %assembla le matrici A, B ed il vettore b
M = buildMassLumping(Me);
Mminus1 = spdiags( 1 ./ diag(M), 0, size(M, 1), size(M, 2));
uu = zeros(size(Me.Nodes.Dof)); 
Dof = Me.Nodes.Dof>0;       %vale 1 per ogni nodo interno

%% Calcolo soluzione stazionaria
uu(Dof) = pcg(A, b, 1e-3, 1000); 
figure;
Me.draw(uu, 'hidemesh');
xlabel('x dir [m]');
ylabel('y dir [m]');
ylabel(colorbar(), 'Displacement [m]');

%% Calcolo delle prime 4 autofunzioni
[V, D] = eigs(Mminus1*A, 4, 'SM');         
%l'ordine crescente o decrescente degli autovalori e dei reEdgesvi 
%autovettori non è assicurato. Devo quindi ordinarli esplicitamente. 
[d, indici] = sort(diag(D));
%indici è un vettore che contiene le permutazioni eseguite tra i vari
%elementi per ordinarli
V=V(:,indici);
figure;
for k = 1:4
    subplot(2, 2, k);
    uu(Dof) = V(:, k); 
    Me.draw(uu, 'hidemesh');
end

%% Analisi dinamica
figure
%Attrio viscoso
a = 1; %attrito viscoso
%Parametri per l'integrazione numerica
dt=0.02;                          %Time step
Tend=5;                          %End time
NumIter=Tend/dt;                  %Number of steps

CentralDisplacement=zeros(NumIter,1);         %vettore di massimi della soluzione nel Time                   
%Inizializzazioni
%valori iniziali di u0 nulli
u0=zeros(size(b));                      

%valori iniziali pari alla prima autofunzione normalizzata a 1
%[V,D]=eigs(A,M,1,'SM');b=b*0;u0=V./max(abs(V(:)));

u1=u0; 
IndexCenter=Me.findClosestNode([0,0]);
%ciclo principale
for k=1:NumIter                           
    
    u2=((2+a*dt)*u1-u0-Mminus1*dt^2*(A*u1-b))/(1+a*dt);
    u0=u1;                      
    u1=u2;
    uu(Dof)=u2;                     %disegno la soluzione
    CentralDisplacement(k)=uu(IndexCenter);
    
    hold off;
    Me.draw(uu,'hidemesh');
    zlim([-1 1]);                         %fisso gli assi in direzione z
    xlabel('x dir [m]');
    ylabel('y dir [m]');
    ylabel(colorbar(),'Displacement [m]');
    title(['t = ' num2str(k*dt) 's']);
    drawnow;
end

figure;
plot((1:NumIter)*dt,CentralDisplacement);
xlabel('t [s]');
ylabel('Displacement [m]');
grid on;
title('Maximum displacement vs. time');
