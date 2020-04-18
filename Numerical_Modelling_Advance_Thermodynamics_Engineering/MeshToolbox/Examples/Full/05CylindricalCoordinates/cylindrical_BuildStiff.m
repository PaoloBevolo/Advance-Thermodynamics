function [D,b] = cylindrical_BuildStiff(Me)
%Assemble the matrix D and the vector b of the Diffusion problem with
%homogeneous B.C.s
%Input:
%   Me     :a Mesh2D object
%
%Output:
%   D      :diffusion matrix
%   b      :constant terms vector

%check inputs
if nargin<2
    f=@(x,y)4;
end 

%for clarity, call some properties of Me with shorter names
V=Me.Triangles.Vertices;
Areas=Me.Triangles.Areas;
CenterOfMass=Me.Triangles.CenterOfMass;
Nodes=Me.Nodes;
Dof=Me.Nodes.Dof;
%number of internal nodes: we know that the N unknown nodes are numbered from
%1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
%(degrees of freedom)
numDof = max(Dof);

%vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
%we save the rows, columns and values corresponding to each contribution; 
%at the end, we'll call sparse(...) to obtain the diffusion matrix
b = zeros(numDof,1);
row = zeros(Me.MatrixContributions,1);
col = zeros(Me.MatrixContributions,1);
d = zeros(Me.MatrixContributions,1);
pos=1;  %we start from the element in position 1, we'll increase this index 
        %everytime we add an entry

mu =Me.evaluateProperty('mu');
beta=Me.evaluateProperty('beta');
%main loop on each triangle        
for e=1:size(V,1)   
    Dx(1) = Nodes.X(V(e,3)) - Nodes.X(V(e,2));
    Dx(2) = Nodes.X(V(e,1)) - Nodes.X(V(e,3));
    Dx(3) = Nodes.X(V(e,2)) - Nodes.X(V(e,1));
    Dr(1) = Nodes.Y(V(e,3)) - Nodes.Y(V(e,2));
    Dr(2) = Nodes.Y(V(e,1)) - Nodes.Y(V(e,3));
    Dr(3) = Nodes.Y(V(e,2)) - Nodes.Y(V(e,1));
    
   
    %we evaluate the external force in the center of mass of this triangle
    rb=CenterOfMass.Y(e);   
    force = f(CenterOfMass.X(e),rb);
    %evaluate the value of the coefficient in front of the Laplace operator

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
                betax=beta(e,1);
                betar=beta(e,2)-mu(e)/rb;
                diffusion=mu(e)*(Dr(ni)*Dr(nj)+Dx(ni)*Dx(nj))/(4.0*Areas(e));
                transport=(betax*Dr(nj)+betar*Dx(nj))*1/6;
                 %%is it unknown as well?
                if jj > 0                                   
                    %Non sparse solution: D(ii,jj)=D(ii,jj) +  diffusion + transport;
                    row(pos)=ii;
                    col(pos)=jj;                    
                    d(pos)= diffusion + transport;
                    pos=pos+1;
                else
                    value=Me.BC.DirichletNodes(-jj,2);
                    b(ii)=b(ii)-(diffusion+transport)*value;
                end
            end
            %build the constant terms vector adding the external
            %contribution
            b(ii) = b(ii) + Areas(e)*force/3.0;
        end
    end
end
%assemble the stiffness matrix D from the
D=sparse(row,col, d, numDof, numDof);