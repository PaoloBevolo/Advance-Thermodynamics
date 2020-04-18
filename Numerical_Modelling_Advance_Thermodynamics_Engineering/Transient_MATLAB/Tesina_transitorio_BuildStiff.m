function [D, bconst, bvar] = heatEquationVariableDirichlet_BuildStiff(Me)
%Assemble the matrix D and the vectors b of the Diffusion problem with 
%time-varying BC.s
%Input:
%   Me     :a Mesh2D object
%
%Output:
%   D      :diffusion matrix
%   b      :constant terms vector

%for clarity, call some properties of Me with shorter names
V=Me.Triangles.Vertices;
Areas=Me.Triangles.Areas*10^3;
Nodes=Me.Nodes;
Dof=Me.Nodes.Dof;
Edges=Me.Edges; %condizione necessario per Robin
Robin=Me.BC.RobinEdges;%condizione necessaria per Robin

%number of internal nodes: we know that the N unknown nodes are numbered from
%1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
%(degrees of freedom)
numDof = max(Dof);

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
bconst = zeros(numDof,1);
bvar = zeros(numDof,1);
row = zeros(Me.MatrixContributions,1);
col = zeros(Me.MatrixContributions,1);
d = zeros(Me.MatrixContributions,1);
pos=1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry

%evaluate the value of the coefficient in front of the Laplace operator
c = Me.evaluateProperty('lambda');
ro_fluid = 975; %densità acqua
cp_fluid = 3950; % cp dell'acqua
beta=((Me.evaluateProperty('beta'))*ro_fluid*cp_fluid); %indica il termine 
        %di trasporto che è aggiunto ma avendo velocità solo nella zona del 
        %tubo posso tenere conto direttamente della densità del fluido del Cp

%main loop on each triangle        
for e=1:size(V,1)   
    Dx(1) = (Nodes.X(V(e,3)) - Nodes.X(V(e,2)))*10^3; %andiamo a scalare la 
                %distanza tra i nodi in modo da considerare l'approssimazione 
                %di 3000m a 3m;
    Dx(2) = (Nodes.X(V(e,1)) - Nodes.X(V(e,3)))*10^3;
    Dx(3) = (Nodes.X(V(e,2)) - Nodes.X(V(e,1)))*10^3;
    Dy(1) = Nodes.Y(V(e,3)) - Nodes.Y(V(e,2));
    Dy(2) = Nodes.Y(V(e,1)) - Nodes.Y(V(e,3));
    Dy(3) = Nodes.Y(V(e,2)) - Nodes.Y(V(e,1));
    
    %for each vertex of this triangle 
    for ni=1:3
        %look at the "unknown" numbering: if the node is positive, it
        %corresponds to a degree of freedom of the problem
        ii = Dof(V(e,ni));
        %is it unknown?
        if ii > 0 
            %yes it is! second loop on the vertices
            for nj=1:3
                jj = Dof(V(e,nj));
                dtmp=c(e)*(Dy(ni)*Dy(nj)+Dx(ni)*Dx(nj))/(4.0*Areas(e)) + ...
                    (-beta(e,1)*Dy(nj)+beta(e,2)*Dx(nj))*1/6; %ho aggiunto 
                                            %il termine relativo a 'beta' 
                %%is it unknown as well?                
                if jj > 0
                    %add the contribution to the stiffness matrix 
                    row(pos)=ii;
                    col(pos)=jj;
                    d(pos)=dtmp;
                    pos=pos+1;
                    %Non sparse solution: D(ii,jj)=D(ii,jj) + c*(Dy(i)*Dy(j)
                    %+Dx(i)*Dx(j))/(4.0*Area) ;
                else %%jj is a node on a Dirichlet edge
                    val=Me.BC.DirichletNodes(-jj,2);              
                    if abs(Me.Nodes.Y(V(e,nj)))>=3.4 %on the inner circle?
                        bvar(ii) = bvar(ii) - dtmp*val ;
                    else
                        bconst(ii) = bconst(ii) - dtmp*val ;
                    end
                end
            end
            %in this example, no external force contribution
        end
    end
end


% Introduco condizione di Robin necessaria per il bordo superiore
for k=1:size(Robin,1)
    Node1=Edges(Robin(k,1),1);
    Node2=Edges(Robin(k,1),2);
    dx=(Nodes.X(Node1)-Nodes.X(Node2))*10^3;
    dy=Nodes.Y(Node1)-Nodes.Y(Node2);    
    dist=sqrt(dx*dx+dy*dy);
    ii1=Dof(Node1);
    ii2=Dof(Node2);
    g=Robin(k,3);
    h=Robin(k,2);
    if ii1>0 && ii2<0 %ii1 is unknown, ii2 is known 
        bconst(ii1)=bconst(ii1)+g/2*dist;
        row(pos)=ii1;
        col(pos)=ii1;
        d(pos)=h*dist/3;
        pos=pos+1;
        %D(ii1,ii1)=D(ii1,ii1)+h*dist/3;
    elseif ii1<0 && ii2>0 %ii1 is known, ii2 is unknown 
        bconst(ii2)=bconst(ii2)+g/2*dist;        
        row(pos)=ii2;
        col(pos)=ii2;
        d(pos)=h*dist/3;
        pos=pos+1;
        %D(ii2,ii2)=D(ii2,ii2)+h*dist/3;
    else  %both are unknwon 
        bconst(ii1)=bconst(ii1)+g/2*dist; %integrale di bordo di Robin
        bconst(ii2)=bconst(ii2)+g/2*dist;
        row(pos:pos+3)=[ii1;ii2;ii1;ii2]; 
        col(pos:pos+3)=[ii1;ii2;ii2;ii1]; 
        d(pos:pos+3)=[2;2;1;1]*h*dist/6; 
        pos=pos+4;
        %D(ii1,ii1)=D(ii1,ii1)+h*dist/3;
        %D(ii2,ii2)=D(ii2,ii2)+h*dist/3;
        %D(ii1,ii2)=D(ii1,ii2)+h*dist/6;
        %D(ii2,ii1)=D(ii2,ii1)+h*dist/6;        
    end
end
%assemble the stiffness matrix D from the
D=sparse(row,col, d, numDof, numDof);