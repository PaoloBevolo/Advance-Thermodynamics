%% TESINA STAZIONARIO GRUPPO 9

%% pulizia del foglio di calcolo
clc
clear all 
close all

%% Discretizzazione dominio: lavoriamo con due sistemi di riferimento diverso 
%tra asse x e asse y (in un caso km e in un caso metri)

S(1) = regions.rectN([0, 0], [3000/1000, 3]); % terreno low
S(2) = regions.rectN([0, 3], [3000/1000, 3.1]); %isolante low
S(3) = regions.rectN([0, 3.1], [3000/1000, 3.4]); %spessore low
S(4) = regions.rectN([0,3.4], [3000/1000, 5.4]); %tubo con raggio interno 1m
S(5) = regions.rectN([0,5.4], [3000/1000, 5.7]); %spessore up
S(6) = regions.rectN([0,5.7], [3000/1000, 5.8]); %isolante up
S(7) = regions.rectN([0, 5.8], [3000/1000, 7.8]); %terreno up
%I comandi S.draw permettono di visualizzare il dominio con i nodi e lati
%numerati.
%figure(); S.draw('e') 
%figure(); S.draw('n')

%% condizioni al contorno

%inseriamo i nodi sul bordo inferiore:
S(1).Borders(1).insertNode(4, [1800/1000 0]);
S(1).Borders(1).insertNode(5, [1200/1000 0]);

%dirichlet bordo inferiore terreno:
S(1).Borders(1).Bc([4 6])=boundaries.dirichlet(12+273);
S(1).Borders(1).Bc(5)=boundaries.dirichlet( @(x,y) ((12-34)/((1200/1000)-...
    (1500/1000))^2)*(x-1500/1000).^2+34+273 ); %dirichlet parabolico 

%dirichlet lato sinistro tubo
S(4).Borders(1).Bc(1) = boundaries.dirichlet(103+273);

%Condizioni Neumann lato destro
for i = 1:length(S)
    S(i).Borders(1).Bc(3) = boundaries.neumann(0);
end

%condizione di Neumann lato sinistro
for i = 1:length(S)
    if i ~= 4
        S(i).Borders(1).Bc(1) = boundaries.neumann(0);
    end
end

% Dati del problema Robin
hc = 5.5; % coeff convettivo per condizione di robin dell'aria
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

figure();S.draw('bc') %Visualiziamo il dominio con le relative condizioni al bordo
title('Dominio e Bc');
ylabel('m');
xlabel('km');

%% Aggiunta delle proprietà
lambda_floor = 2; %valore di lambda del terreno
lambda_fluid = 0.5; %valore di lambda del fluido
lambda_isulator = 0.04; %valore di lambda dell'isolante
lambda_metal = 50; %valore di lambda del tubo metallico

ro_fluid = 975; %densità acqua
cp_fluid = 3950; % cp dell'acqua

S(4).addProperty('lambda',lambda_fluid); %diffusibilità termica su dominio del tubo
S(1).addProperty('lambda',lambda_floor); %proprietà di diffusione del terreno
S(7).addProperty('lambda',lambda_floor); %proprietà di diffusione del terreno
S(2).addProperty('lambda',lambda_isulator); %proprietà di diffusione dell'isolante
S(6).addProperty('lambda',lambda_isulator); %proprietà di diffusione dell'isolante
S(3).addProperty('lambda',lambda_metal); %proprietà di diffusione del tubo
S(5).addProperty('lambda',lambda_metal); %proprietà di diffusione del tubo

%% Mesh
Me_dimention = [0.1 0.01 0.005]; % Diverse dimensioni dell'area massima della mesh
zz = zeros(100,length(Me_dimention)); %Inizializzazione vettore per calcolo 
                                                            %convergenza mesh

for i=1:length(Me_dimention) %Ciclo per ciascuna dimensione dell'area massima 
                                                                %della mesh
    
Me_S = mesh2D(S,Me_dimention(i)); %Rappresentiamo la mesh
figure();  Me_S.draw();
title(['Mesh area max : ',num2str(Me_dimention(i)),]);
ylabel('m');
xlabel('km');

%% Calcolo della velocità: oltre a valutarne il valore si è deciso di 
%rappresentare il suo andamento 
Me_Tubo = Me_S.extractMesh(4); % estraggo la mesh che mi serve ossia la parte 4
f_beta = @(y) -(0.5)*((y-4.4).^2)+(0.5); % funzione parabolica della 
                                            %velocità assegnata da testo
y_b = Me_Tubo.Nodes.Y; % si calcolano i nodi lungo l'asse y
b_x = f_beta(y_b) ; %si definisce la velocità lungo asse x
b_y = zeros(length(y_b),1); %si definisce la velocità lungo asse y che è 
                        %nulla ed è definita come un vettore colonna di zeri
figure(); 
axis equal; 
quiver(Me_Tubo.Nodes.X,Me_Tubo.Nodes.Y,b_x,b_y); % si disegna il profilo delle velocità 
title(['Profilo di velocità con mesh di area max : ',num2str(Me_dimention(i)),]);
ylabel('m');
xlabel('km');
%Definizione del vettore velocità nelle sue due componenti
XX=Me_Tubo.Triangles.CenterOfMass.X; 
YY=Me_Tubo.Triangles.CenterOfMass.Y; 
beta_x=(Me_Tubo.interpolate(b_x,[XX,YY],1:length(XX)))*10^-3; % si calcola 
        %la componente x della velocità per ciuascun centro di massan 
        %( il comando 1:length(XX) mi permette di valutare tutti i centri
        %di massa della zona indicata)
beta_y=Me_Tubo.interpolate(b_y,[XX,YY],1:length(YY)); %si calcola la componente 
                                                        %y della velocità
V=[beta_x , beta_y]; %si definisce il vettore velocità
%Assegnazione della propriètà di velocità nel condotto
S.addProperty('beta',[0 0]); % si definisce su tutto il dominio beta nullo
S(4).addProperty('beta',[V(:,1) V(:,2)]); %si sovrascive la velocità nel tubo


%% risoluzione del sistema con function 
[A,b] = Tesina_stazionario_BuildStiff(Me_S, @(x,y) zeros(size(x))); %si impone forzante nulla
% Risoluzione del sistema
t = A\b;
T = Me_S.copyToAllNodes(t); %si restituisce il vettore t esteso a tutti i 
                            %nodi e non indicizzato rispetto ai soli gradi 
                            %di libertà

%% Rappresentazione grafica
figure();
axis equal; 
draw(Me_S,T,'h');
ylabel('m');
xlabel('km');
zlabel('K');
title(['Mesh area max :',num2str(Me_dimention(i))]);

%% Elaborazione della soluzione per calcolo convergenza
yy = 4.4*ones(100,1);
xx = linspace(0,3,100).';
zz(:,i) = Me_S.interpolate(T,[xx,yy]);

end

%%Si plotta la convergenza
figure();
axis equal;
plot(xx,zz(:,1:end));
ylabel('m');
xlabel('km');
legend([num2str(Me_dimention(1,1))],[num2str(Me_dimention(1,2))],[num2str(Me_dimention(1,3))]);
title(['Temperatura: convergenza della mesh']);

%% Si valuta il numero di Peclet (non espressamente richiesto dai dati)

Areas=(Me_S.Triangles.Areas)*(10^3); %Si inserisce il vettore che continene 
                                    %il numero delle aree di ciascun triangolo
Num_peclet = zeros(length(Areas),1); %si inizializza Peclet a zeri. Questo 
                                    %vettore deve essere grande nello stesso 
                                    %numero dei triangoli che compongono la mesh
lato = zeros(length(Areas),1); %contiene la lunghezza del lato di ciascun 
                                %triangolo che compone la mesh e lo si 
                                %inizializza a zero
Beta = Me_S.evaluateProperty('beta'); %contiene la proprietà relativa alla velocità
Lambda = Me_S.evaluateProperty('lambda'); %contiene l'informazione legata alla conducibilità
for k = 1:length(Areas) %ciclo triangolo per triangolo
    lato(k) = ((4.*Areas(k))/(sqrt(3))).^0.5; %tramite la formula inversa 
                                    %dell'area mi calcolo la lunghezza del 
                                    %lato con l'ipotesi che i triangoli siano
                                    %tutti equilateri
    Num_peclet(k) = (Beta(k,1)*ro_fluid*cp_fluid*lato(k))/(6*Lambda(k)); %calcolo
                                                %peclet per ciascun triangolo
end
Peclet_max = max(Num_peclet) %Peclet rappresenta il numero massimo dei peclet 
                                %che ho nella triangolazione


%% Preconditioning study
%figure; spy(A);
[L,U] = lu(A);
%figure; spy(L);
%figure; spy(U);
% Inizializzazione dei parametri per il ciclo
opts.type = 'ilutp'; % Threshold dropping e pivoting
Nrep = 10; % Numero di ripetizioni su cui mediare il tempo
droptols = 10.^(-2:-0.5:-8); % Vettore delle tolleranze considerate
% Allocazione dei vettori
time_lu = zeros(size(droptols));
time_bicg = zeros(size(droptols));
flag = zeros(size(droptols));
res = zeros(size(droptols));
iter = zeros(size(droptols));
for k=1:length(droptols)
    opts.droptol = droptols(k);
    time_lu(k) = 0;
    time_bicg(k) = 0;
    for i=1:Nrep
        % Fattorizzazione LU incompleta
        tic();
        [L,U] = ilu(A,opts);
        time_lu(k) = time_lu(k) + toc();
        % Gradiente BiConiugato precondizionato
        tic();
        [x,flag(k),res(k),iter(k)] = bicg(A,b,1e-6,100,L,U);
        time_bicg(k) = time_bicg(k) + toc();
    end
    time_lu(k) = time_lu(k)/Nrep;
    time_bicg(k) = time_bicg(k)/Nrep;
end
% Rappresentazione grafica dei tempi parziali e del tempo totale
figure();
semilogx(droptols,time_lu, droptols,time_bicg, droptols,time_lu+time_bicg);
legend('ilu','bicg','totale');
ylabel('time [s]');
xlabel('droptolerance');
title(['Preconditioning study']);
