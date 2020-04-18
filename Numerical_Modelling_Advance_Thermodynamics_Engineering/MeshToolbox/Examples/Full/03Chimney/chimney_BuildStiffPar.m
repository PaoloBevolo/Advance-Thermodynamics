function [D, bconst, bvar] = chimney_BuildStiff(Me)
%Assemble the matrix D and the vector b of the Diffusion problem with
%Robin B.C.s and non homogeneous Dirichlet B.C.s
%Input:
%   Me     :a Mesh2D object
%   f      :MATLAB function of (x,y) which returns the values of the
%           external source. Default: constant value=4
%
%Output:
%   D      :diffusion matrix
%   bconst :constant terms vector, time independent components
%   bvar   :constant terms vector, time dependent components

%for clarity, call some properties of Me with shorter names
V=Me.Triangles.Vertices;
Dof=Me.Nodes.Dof;
Nodes=Me.Nodes;
Areas=Me.Triangles.Areas;

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


%main loop on each triangle  

    V=Me.Triangles.Vertices;
    Dof=Me.Nodes.Dof;
    Nodes=Me.Nodes;
    Areas=Me.Triangles.Areas;

    for e=1:size(V,1)   

    Dx=[
    Nodes.X(V(e,3)) - Nodes.X(V(e,2));
    Nodes.X(V(e,1)) - Nodes.X(V(e,3));
    Nodes.X(V(e,2)) - Nodes.X(V(e,1))];
    Dy=[
    Nodes.Y(V(e,3)) - Nodes.Y(V(e,2));
    Nodes.Y(V(e,1)) - Nodes.Y(V(e,3));
    Nodes.Y(V(e,2)) - Nodes.Y(V(e,1))];
    
    Area=Areas(e);
    %evaluate the value of the coefficient in front of the Laplace operator
    c=Me.evaluateProperty('c',e);
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
                dtmp=c*(Dy(ni)*Dy(nj)+Dx(ni)*Dx(nj))/(4.0*Area) ;
                %is it unknown as well?
                if jj > 0
                    %add the contribution to the stiffness matrix 
                    tmp=pos;
                    row(tmp)=ii;
                    col(tmp)=jj;                    
                    d(tmp)=dtmp;
                    pos=pos+1;
                else 
                    %Non homogeneous Dirichlet B.C.s
                    val=Me.BC.DirichletNodes(-jj,2);     
                    xPositionNodeJ=Nodes.X(V(e,nj),1);
                    if xPositionNodeJ==-0.5
                        bconst(ii) = bconst(ii) - dtmp*val ;
                    else
                        bvar(ii) = bvar(ii) - dtmp*val ;                        
                    end
                end
            end
        end
    end
    end

end

Dof=Me.Nodes.Dof;
Nodes=Me.Nodes;

%%%%%%%%%%%%% ROBIN B.C.
Edges=Me.Edges;
Robin=Me.BC.RobinEdges;

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
        bconst(ii1)=bconst(ii1)+g/2*dist;
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
%assemble the stiffness matrix D
D=sparse(row, col, d, numDof,numDof);
