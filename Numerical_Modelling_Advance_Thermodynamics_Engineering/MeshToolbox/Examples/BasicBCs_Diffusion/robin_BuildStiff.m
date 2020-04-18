function [D,b] = robin_BuildStiff(Me,f)
%Assemble the matrix D and the vector b of the Diffusion problem with
%Robin B.C.s
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
Robin=Me.BC.RobinEdges;
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
                %is it unknown as well?
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

for k=1:size(Robin,1)
    Node1=Edges(Robin(k,1),1);
    Node2=Edges(Robin(k,1),2);
    dx=Nodes.X(Node1)-Nodes.X(Node2);
    dy=Nodes.Y(Node1)-Nodes.Y(Node2);    
    dist=sqrt(dx*dx+dy*dy);
    ii1=Dof(Node1);
    ii2=Dof(Node2);
    g=Robin(k,3);
    h=Robin(k,2);
    if ii1>0 && ii2<0 %ii1 is unknown, ii2 is known
        b(ii1)=b(ii1)+g/2*dist;
        row(pos)=ii1;
        col(pos)=ii1;
        d(pos)=h*dist/3;
        pos=pos+1;
        %D(ii1,ii1)=D(ii1,ii1)+h*dist/3;
    elseif ii1<0 && ii2>0 %ii1 is known, ii2 is unknown
        b(ii2)=b(ii2)+g/2*dist;        
        row(pos)=ii2;
        col(pos)=ii2;
        d(pos)=h*dist/3;
        pos=pos+1;
        %D(ii2,ii2)=D(ii2,ii2)+h*dist/3;
    else  %both are unknwon
        b(ii1)=b(ii1)+g/2*dist;
        b(ii2)=b(ii2)+g/2*dist;
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
%assemble the stiffness matrix D
D=sparse(row, col, d, numDof, numDof);