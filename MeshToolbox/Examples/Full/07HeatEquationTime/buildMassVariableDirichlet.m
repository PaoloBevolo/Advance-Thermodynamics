function [M, mvar] = buildMassVariableDirichlet(Me)
%Assemble the mass matrix M
%Input:
%   Me     :a Mesh2D object
%
%Output:
%   M      :mass matrix
%   mvar   :vector to include the effect of the non Homogeneous Dirichlet B.C. 

%for clarity, call some properties of Me with shorter names
V = Me.Triangles.Vertices;
Dof = Me.Nodes.Dof;
Areas=Me.Triangles.Areas;
%number of internal nodes: we know that the N unknown nodes are numbered from
%1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
%(degrees of freedom)
numDof = max(Dof);

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
row = zeros(Me.MatrixContributions, 1);
col = zeros(Me.MatrixContributions, 1);
m = zeros(Me.MatrixContributions, 1);
mvar = zeros(numDof, 1);
pos = 1; %we start from the element in position 1, we'll increase this index 
       %everytime we add an entry

%evaluate the value of the coefficient in front of the time derivative
rho = Me.evaluateProperty('rho');
%main loop on each triangle  
for e = 1:size(V, 1)
    %for each vertex of this triangle 
    for ni = 1:3
        %look at the "unknown" numbering: if the node is positive, it
        %corresponds to a degree of freedom of the problem
        ii = Dof(V(e, ni));
         %is it unknown?
        if ii > 0    
            %yes it is! second loop on the vertices
            for nj = 1:3
                jj = Dof(V(e, nj));
                %%is it unknown as well?
                mtmp = 1/12 * Areas(e) * rho(e) * ((ii == jj) + 1);
                if jj > 0
                    %add the contribution to the mass matrix                     
                    row(pos) = ii;
                    col(pos) = jj;
                    m(pos) = mtmp;
                    pos = pos + 1;
                else
                    if abs(Me.Nodes.X(V(e,nj)))<=0.2 %on the inner circle?
                        mvar(ii)=mvar(ii)+mtmp;
                    end
                end
            end
        end
    end
end
%assemble the mass matrix M
M = sparse(row, col, m, numDof, numDof);
