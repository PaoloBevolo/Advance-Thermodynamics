function [D,R,b] = coupledDirichlet_BuildStiff(Me)
%Builds the stiffness matrices for a diffusion-reaction problem
%Input:
%   Me    :a Mesh2D object
%
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
%Number of unknown nodes
numDof = max(Dof);

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
b = zeros(numDof,1);
row = zeros(Me.MatrixContributions,1);
col = zeros(Me.MatrixContributions,1);
d = zeros(Me.MatrixContributions,1);
r = zeros(Me.MatrixContributions,1);
pos=1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry
%evaluate the value of the coefficient in front of the Laplace operator
mu=Me.evaluateProperty('mu');
%evaluate the value of the couling coefficient 
sigma=Me.evaluateProperty('sigma');%main loop on each triangle        
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
        ii = Dof(V(e,ni));
        %is it unknown?
        if ii > 0 
            %yes it is! second loop on the vertices
            for nj=1:3
                jj = Dof(V(e,nj));                
                dtmp=mu(e)*(Dy(ni)*Dy(nj)+Dx(ni)*Dx(nj))/(4.0*Areas(e));
                rtmp=sigma(e)*Areas(e)*(1+(ni==nj))/12;
                %is it unknown as well?
                if jj > 0
                     %add the contribution to the stiffness & reaction matrices 
                    row(pos)=ii;
                    col(pos)=jj;
                    d(pos)=dtmp;
                    r(pos)=rtmp;
                    pos=pos+1;
                else
                   %build the constant terms vector 
                    val=Me.BC.DirichletNodes(-jj,2);
                    b(ii) = b(ii) - (dtmp+rtmp)*val;
                end
            end
        end
    end
end
%assemble the stiffness matrix D and the reaction matrix R
D=sparse(row,col, d,numDof,numDof);
R=sparse(row,col, r,numDof,numDof);

