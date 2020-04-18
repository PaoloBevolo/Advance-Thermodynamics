function [D,b] = neumannNonHomo_BuildStiff(Me,f)
%Assemble the matrix D and the vector b of the Diffusion problem with
%Neumann non homogeneous B.C.s
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
    f=@(x,y)zeros(size(x))+4;
end 

%for clarity, call some properties of Me with shorter names
V=Me.Triangles.Vertices;
Areas=Me.Triangles.Areas;
CenterOfMass=Me.Triangles.CenterOfMass;
Nodes=Me.Nodes;
Dof=Me.Nodes.Dof;
Edges=Me.Edges;
Neumann=Me.BC.NeumannEdges;
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

%we evaluate the external force in the center of mass of this triangle
force = f(CenterOfMass.X,CenterOfMass.Y);
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
                %%is it unknown as well?
                if jj > 0
                    %add the contribution to the stiffness matrix 
                    row(pos)=ii;
                    col(pos)=jj;
                    d(pos)=c(e)*(Dy(ni)*Dy(nj)+Dx(ni)*Dx(nj))/(4.0*Areas(e));
                    pos=pos+1;
                    %Non sparse solution: D(ii,jj)=D(ii,jj) + c*(Dy(i)*Dy(j)+Dx(i)*Dx(j))/(4.0*Area) ;
                end
            end
            %build the constant terms vector adding the external
            %contribution
            b(ii) = b(ii) + Areas(e)*force(e)/3.0;
        end
    end
end

for k=1:size(Neumann,1)
    Node1=Edges(Neumann(k),1);
    Node2=Edges(Neumann(k),2);
    dx=Nodes.X(Node1)-Nodes.X(Node2);
    dy=Nodes.Y(Node1)-Nodes.Y(Node2);    
    dist=sqrt(dx*dx+dy*dy);
    g=Me.BC.NeumannEdges(k,2);
    if Dof(Node1)>0 
        b(Dof(Node1))=b(Dof(Node1))+g/2*dist;
    end
    if Dof(Node2)>0 
        b(Dof(Node2))=b(Dof(Node2))+g/2*dist;
    end
end
%assemble the stiffness matrix D from the
D=sparse(row, col, d, numDof, numDof);