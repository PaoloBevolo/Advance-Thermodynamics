function [D,b] = dirichletHomo_BuildStiff(Me,f)
%Assemble the matrix D and the vector b of the Diffusion problem with
%homogeneous B.C.s
%Input:
%   Me     :a Mesh2D object
%   f      :MATLAB function of (x,y) which returns the values of the
%           external source. Default: constant value=4
%
%Output:
%   D      :diffusion matrix
%   b      :constant terms vector

%check inputs
if nargin<2
    f = @(x,y)zeros(size(x))-4;
end 

%for clarity, call some properties of Me with shorter names
%Indices of vertices of each triangle; 3 cols. matrix
V = Me.Triangles.Vertices;		
%Area of each triangle; column vector
Areas = Me.Triangles.Areas;		
%x and y coordinates of the centres of mass;col. vect.
CenterOfMass = Me.Triangles.CenterOfMass; 
%x and y coordinates of the certices; column vectors
X=Me.Nodes.X; Y=Me.Nodes.Y; 
%Numbering of Dirichlet(<0) and unknown (>0) nodes
Dof = Me.Nodes.Dof;
%number of internal nodes: we know that the numDof unknown nodes are 
%numbered from 1 to numDof in Me.Nodes.Dof; its maximum is therefore the 
%number of unknown (degrees of freedom)
numDof = max(Dof);
 
%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix

row = zeros(Me.MatrixContributions, 1);
col = zeros(Me.MatrixContributions, 1);
d = zeros(Me.MatrixContributions, 1);
pos = 1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry
b = zeros(numDof, 1);
      
%external force, evaluated in the center of mass of each triangle
force = f(CenterOfMass.X,CenterOfMass.Y);

%evaluate the value of the coefficient in front of the Laplace operator
mu = Me.evaluateProperty('mu');

%main loop on each triangle        
for e= 1:size(V, 1)   
    Dx(1) = X(V(e, 3)) - X(V(e, 2));
    Dx(2) = X(V(e, 1)) - X(V(e, 3));
    Dx(3) = X(V(e, 2)) - X(V(e, 1));
    Dy(1) = Y(V(e, 3)) - Y(V(e, 2));
    Dy(2) = Y(V(e, 1)) - Y(V(e, 3));
    Dy(3) = Y(V(e, 2)) - Y(V(e, 1));
    
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
                if jj > 0
                    %add the contribution to the stiffness matrix 
                    row(pos) =ii;
                    col(pos) = jj;
                    d(pos) = mu(e) * (Dy(ni) * Dy(nj) + Dx(ni) * Dx(nj)) / (4 * Areas(e));
                    pos = pos + 1;
                    %Non sparse solution: D(ii,jj)=D(ii,jj) + c*(Dy(i)*Dy(j)+Dx(i)*Dx(j))/(4.0*Area) ;
                end
            end
            %build the constant terms vector adding the external contribution            
            b(ii) = b(ii) + Areas(e) * force(e) / 3;
        end
    end
end
%assemble the stiffness matrix D
D = sparse(row, col, d, numDof, numDof);