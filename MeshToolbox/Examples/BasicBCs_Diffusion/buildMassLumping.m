function M = buildMassLumping(Me)
%Assemble the diagonal mass lumped matrix M, stored as a vector
%Input:
%   Me     :a Mesh2D object
%
%Output:
%   M      :sparse diagonal mass matrix

%for clarity, call some properties of Me with shorter names
V = Me.Triangles.Vertices;
Dof = Me.Nodes.Dof;
Areas=Me.Triangles.Areas;
%number of internal nodes: we know that the N unknown nodes are numbered from
%1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
%(degrees of freedom)
numDof = max(Dof);

%vectors preallocation: I store the data of the main diagonal in  a vector
m = zeros(numDof,1);
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
            m(ii) = m(ii) + Areas(e) * rho(e) / 3;
        end
    end
end
%assemble the diagonal mass matrix M
M = spdiags(m, 0, numDof, numDof);