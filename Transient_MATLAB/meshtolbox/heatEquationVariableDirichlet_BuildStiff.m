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
Areas=Me.Triangles.Areas;
Nodes=Me.Nodes;
Dof=Me.Nodes.Dof;
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
c = Me.evaluateProperty('mu');

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
        ii = Dof(V(e,ni));
        %is it unknown?
        if ii > 0 
            %yes it is! second loop on the vertices
            for nj=1:3
                jj = Dof(V(e,nj));
                dtmp=c(e)*(Dy(ni)*Dy(nj)+Dx(ni)*Dx(nj))/(4.0*Areas(e));
                %%is it unknown as well?                
                if jj > 0
                    %add the contribution to the stiffness matrix 
                    row(pos)=ii;
                    col(pos)=jj;
                    d(pos)=dtmp;
                    pos=pos+1;
                    %Non sparse solution: D(ii,jj)=D(ii,jj) + c*(Dy(i)*Dy(j)+Dx(i)*Dx(j))/(4.0*Area) ;
                else %%jj is a node on a Dirichlet edge
                    val=Me.BC.DirichletNodes(-jj,2);              
                    if abs(Me.Nodes.X(V(e,nj)))<=0.2 %on the inner circle?
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
%assemble the stiffness matrix D from the
D=sparse(row,col, d, numDof, numDof);