%% TESINA GRUPPO 9 TRANSIENTE
%% pulizia rimanenze
clc
clear all 
close all

%% Discretizzazione dominio: lavoriamo con due sistemi di riferimento 
            %diverso tra asse x e asse y (in un caso km e in un caso metri)

S(1) = regions.rectN([0, 0], [3000/1000, 3]); % terreno low
S(2) = regions.rectN([0, 3], [3000/1000, 3.1]); %isolante low
S(3) = regions.rectN([0, 3.1], [3000/1000, 3.4]); %spessore low
S(4) = regions.rectN([0,3.4], [3000/1000, 5.4]); %tubo con raggio interno 1m
S(5) = regions.rectN([0,5.4], [3000/1000, 5.7]); %spessore up
S(6) = regions.rectN([0,5.7], [3000/1000, 5.8]); %isolante up
S(7) = regions.rectN([0, 5.8], [3000/1000, 7.8]); %terreno up

%% condizioni al contorno

%inseriamo i nodi sul bordo inferiore:
S(1).Borders(1).insertNode(4, [1800/1000 0]);
S(1).Borders(1).insertNode(5, [1200/1000 0]);

%dirichlet bordo inferiore terreno:
S(1).Borders(1).Bc([4 6])=boundaries.dirichlet(12+273);
S(1).Borders(1).Bc(5)=boundaries.dirichlet( @(x,y) ((12-34)/((1200/1000)-...
    (1500/1000))^2)*(x-1500/1000).^2+34+273 ); %dirichlet parabolico 

%dirichlet lato sinistro tubo
S(4).Borders(1).Bc(1) = boundaries.dirichlet(1); %impongo 1 per poter 
                                       %andare a inserire l'andamento 
                                       %di temperatura in funzione del tempo

%Condizioni Neumann lato destro
for i = 1:length(S)
    S(i).Borders(1).Bc(3) = boundaries.neumann(0);
end

%condizione di Neuman lato sinistro
for i = 1:length(S)
    if i ~= 4
        S(i).Borders(1).Bc(1) = boundaries.neumann(0);
    end
end

% Dati del problema robin
hc = 5.5; %per coeff convettivo per condizione di robin dell'aria
t_inf = 20+273; % T infinito di robin

%condizione di Robin
a = hc;
b = hc*t_inf;
S(7).Borders(1).Bc(2)=boundaries.robin(a,b);

%condizione di continuità
for i = 1:length(S)-1
    S(i).Borders(1).Bc(2)=boundaries.continuity();
end

for i = 2:length(S)
    S(i).Borders(1).Bc(4)=boundaries.continuity();
end

%% Aggiunta delle proprietà
lambda_floor = 2; %valore di lambda del terreno
lambda_fluid = 0.5; %valore di lambda del fluido
lambda_isulator = 0.04; %valore di lambda dell'isolante
lambda_metal = 50; %valore di lambda del tubo metallico

cp_fluid = 3950; % cp dell'acqua
cp_metal = 500; %cp dell'acciaio
cp_isulator = 1300; %cp dell'isolante 
cp_floor = 900; %cp del terreno

ro_fluid = 975; %densità acqua
ro_metal = 7800; %densità metallo
ro_isulator = 250; %densità isolante
ro_floor = 1100; %densità del terreno

S(4).addProperty('lambda',lambda_fluid); %diffusibilità termica su dominio del tubo
S(1).addProperty('lambda',lambda_floor); %proprietà di diffusione del terreno
S(7).addProperty('lambda',lambda_floor); %proprietà di diffusione del terreno
S(2).addProperty('lambda',lambda_isulator); %proprietà di diffusione dell'isolante
S(6).addProperty('lambda',lambda_isulator); %proprietà di diffusione dell'isolante
S(3).addProperty('lambda',lambda_metal); %proprietà di diffusione del tubo
S(5).addProperty('lambda',lambda_metal); %proprietà di diffusione del tubo

S(4).addProperty('rho',ro_fluid*cp_fluid); %cp * ro su dominio del tubo
S(1).addProperty('rho',ro_floor*cp_floor); %cp * ro del terreno
S(7).addProperty('rho',ro_floor*cp_floor); %cp * ro del terreno
S(2).addProperty('rho',ro_isulator*cp_isulator); %cp * ro dell'isolante
S(6).addProperty('rho',ro_isulator*cp_isulator); %cp * ro dell'isolante
S(3).addProperty('rho',ro_metal*cp_metal); %cp * ro del tubo
S(5).addProperty('rho',ro_metal*cp_metal); %cp * ro del tubo


%% Mesh
Me_dimention = [0.005]; %valore della dimensione della meh 
Me_S = mesh2D(S,Me_dimention);

%% Calcolo della velocità: 
Me_Tubo = Me_S.extractMesh(4); % estraggo la mesh che mi serve ossia la parte 4
f_beta = @(y) -(0.5)*((y-4.4).^2)+(0.5); % funzione parabolica della velocità
y_b = Me_Tubo.Nodes.Y; % calcolo i nodi lungo l'asse y
b_x = f_beta(y_b) ; %definisco la velocità lungo asse x
b_y = zeros(length(y_b),1); %definisco velocità lungo asse y che è nula ed è 
                                %definita come un vettore colonna di zeri

XX=Me_Tubo.Triangles.CenterOfMass.X; % calcolo i centri di massa di x per 
                                    %ciarcun triangolo per la cordinata x
YY=Me_Tubo.Triangles.CenterOfMass.Y; % calcolo i centri di massa di x per 
                                    %ciarcun triangolo per la cordinata x
beta_x=(Me_Tubo.interpolate(b_x,[XX,YY],1:length(XX)))*10^-3; % calcolo la 
                        %componente x della velocità per ciuascun centro di 
                        %massa ( il comando 1:length(XX) mi permette di 
                        %valutare tutti i centri di massa della zona indicata)
beta_y=Me_Tubo.interpolate(b_y,[XX,YY],1:length(YY)); %calcolo la componente y 
                                                            %della velocità
V=[beta_x , beta_y]; %definisco il vettore velocità

%% definisco la proprietà di velocità:

S.addProperty('beta',[0 0]); % definisco su tutto il dominio beta nullo
S(4).addProperty('beta',[V(:,1) V(:,2)]); %sovrascrivo velocità diversa da 
                                %zero nel dominio del tubo tale per cui ho 
                                %V(:,1)=tutti i valori della prima coordinata 
                                %di V e V(:,2)= tutti i valori della seconda 
                                %coordinata di V

%% risoluzione del sistema con function 
[A,bconst_iniziale, bvar_iniziale] = Tesina_transitorio_INIZIALE_BuildStiff(Me_S); 
                        %richiamo la funzione per il calcolo della condizione 
                        %iniziale

[D,bconst, bvar] = Tesina_transitorio_BuildStiff(Me_S);%richiamo la funzione per 
                        %il calcolo di tutti gli altri istanti temporali

%% Nodi sul bordo di dirichlet annullando il profilo parabolico
DirichletNodesVar=Me_S.find(@(x,y) y>=3.4,'d'); %definisco i nodi sul lato 
                                                    %di dirichlet che varia
%% Gradi di libertà
Dof=Me_S.Nodes.Dof>0; %definisco i gradi di libertà
%% Soluzione stazionaria
T0=12+273; %temperatura iniziale
uStationary0=A\(bconst_iniziale+bvar_iniziale*T0);
uu=Me_S.copyToAllNodes(uStationary0);
uu(DirichletNodesVar)=T0;
figure; 
Me_S.draw( uu);
ylabel(colorbar(),'Temperature [K]');

%% time evolution: definition of the function T
%edge 4 temperature
Tend=20000000;
fDirichlet=@(t) 273+12+(91)*(2/pi)*atan(t);
[M, mvar] = Tesina_Massa(Me_S);

%% implicit Euler method
disp('Implicit Euler method');
u=uStationary0;
dt=25000;
figure;
S=(M+D*dt);
tic
TPipeOld=fDirichlet(0);

TT=zeros(size(Tend));
tt=zeros(size(Tend));

for k=1:Tend/dt, 
    t=k*dt;
    TPipe=fDirichlet(t);
    %original: u=(M+D*dt)\(M*u+b*dt*T4);
    u=S\(M*u+dt*(bconst+bvar*TPipe)-mvar*(TPipe-TPipeOld));
    uu(Dof)=u;
    uu(DirichletNodesVar)=TPipe;
    TPipeOld=TPipe;
    hold off;  
    TT(k) = uu(1020);
    tt(k) = t;
    if rem (k,100)==0
        Me_S.draw(uu,'hidemesh');
        zlim([273 380]);
        caxis([273 376]);
        title(['Andamento termico al tempo t = ' num2str(t) 's']);    
        xlabel('x dir [m]');
        ylabel('y dir [km]');
        zlabel('[K]');
        ylabel(colorbar(),'Temperature [K]');
        drawnow();
    end
end 
toc

figure();
plot(tt, TT);
title('andamento temperatura nodo 1020');
xlabel('time [s]');
ylabel('Temperature [K]');

%% Rappresentazione dei primi 2000000 secondi
Tend=2000000;
u=uStationary0;
dt=1000;
figure;
S=(M+D*dt);
tic
TPipeOld=fDirichlet(0);

for k=1:Tend/dt, 
    t=k*dt;
    TPipe=fDirichlet(t);
    %original: u=(M+D*dt)\(M*u+b*dt*T4);
    u=S\(M*u+dt*(bconst+bvar*TPipe)-mvar*(TPipe-TPipeOld));
    uu(Dof)=u;
    uu(DirichletNodesVar)=TPipe;
    TPipeOld=TPipe;
    hold off;  
    if rem (k,100)==0
        Me_S.draw(uu,'hidemesh');
        zlim([273 380]);
        caxis([273 376]);
        title(['Dettaglio iniziale della temperatura per t = ' num2str(t) 's']);    
        xlabel('x dir [m]');
        ylabel('y dir [km]');
        zlabel('[K]');        
        ylabel(colorbar(),'Temperature [K]');
        drawnow();
    end
end
toc

