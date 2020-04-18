function [D,R,b] = coupledNeumann_BuildStiff(Me)
%coupledNeumann_BuildStiff builds the stiffness matrices for a diffusion-reaction problem
%The corresponding linear system will return the solution for ALL the
%nodes, not only the unknown ones. 
%Input:
%   Me   :  mesh2d object
%shape properties:
%   c, a
%Output:
%   D      :stiffness matrix - diffusion terms
%   R      :stiffness matrix - reaction terms
%   b      :constant terms vector
%
%Copyright 2014 Paolo Bardella

%For clarity, I call some variables with shorter names
V=Me.Triangles.Vertices;
Dof=Me.Nodes.Dof;
Nodes=Me.Nodes;
Areas=Me.Triangles.Areas;
%Total number of nodes
NumTotalNodes = length(Dof);

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
NumDirichletNodes=size(Me.BC.DirichletNodes,1);
row=zeros(Me.MatrixContributions+NumDirichletNodes,1);
col=zeros(Me.MatrixContributions+NumDirichletNodes,1);
r = zeros(Me.MatrixContributions+NumDirichletNodes,1);
d = zeros(Me.MatrixContributions+NumDirichletNodes,1);
b = zeros(NumTotalNodes,1);
pos=1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry

%evaluate the value of the coefficient in front of the Laplace operator
c = Me.evaluateProperty('mu');
%evaluate the value of the couling coefficient
a=Me.evaluateProperty('sigma');
%main loop on each triangle        
for e=1:size(V,1)   
    Dx(1) = Nodes.X(V(e,3)) - Nodes.X(V(e,2));
    Dx(2) = Nodes.X(V(e,1)) - Nodes.X(V(e,3));
    Dx(3) = Nodes.X(V(e,2)) - Nodes.X(V(e,1));
    Dy(1) = Nodes.Y(V(e,3)) - Nodes.Y(V(e,2));
    Dy(2) = Nodes.Y(V(e,1)) - Nodes.Y(V(e,3));
    Dy(3) = Nodes.Y(V(e,2)) - Nodes.Y(V(e,1));


    
    %for each vertex of this triangle 
    for ni=1:3
        %look at the "unknown" numbering: if the node is positive, it
        %corresponds to a degree of freedom of the problem
        ii = V(e,ni);
        %is it unknown?
        if Dof(ii) > 0
            %yes it is! second loop on the vertices
            for nj=1:3
                jj = V(e,nj);                
                dtmp=c(e)*(Dy(ni)*Dy(nj)+Dx(ni)*Dx(nj))/(4.0*Areas(e));
                rtmp=a(e)*Areas(e)*(1+(ni==nj))/12;   
                %%is it unknown as well?
                if Dof(jj) > 0
                    row(pos)=ii;
                    col(pos)=jj;
                    d(pos)=dtmp;
                    r(pos)=rtmp;
                    pos=pos+1;                    
                else
                    b(ii)=b(ii)-(dtmp+rtmp)*Me.BC.DirichletNodes(-Dof(jj),2);
                end
                
            end
        end
    end
end

%Contributions from Dirichlet nodes: 1s on the main  diagonal of D or R (here I choose D) 
for k=1:NumDirichletNodes
    ii=Me.BC.DirichletNodes(k,1);
    row(pos)=ii;
    col(pos)=ii;
    d(pos)=1;
    r(pos)=0;    
    pos=pos+1;
    b(ii)=Me.BC.DirichletNodes(k,2);
end

%assemble the stiffness matrix D and the reaction matrix R
D=sparse(row, col, d, NumTotalNodes,NumTotalNodes);
R=sparse(row, col, r, NumTotalNodes,NumTotalNodes);