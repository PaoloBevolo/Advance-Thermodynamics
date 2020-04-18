function [Dfinal,bfinal] = dirichletHomo_BuildStiffPar(Me,f, NPools)
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
    f = @(x,y)4;
end 
if nargin<3
    Pool=gcp;
    NPools=Pool.NumWorkers;
end
    

%for clarity, call some properties of Me with shorter names
V = Me.Triangles.Vertices;
Areas = Me.Triangles.Areas;
CenterOfMass = Me.Triangles.CenterOfMass;
Nodes = Me.Nodes;
Dof = Me.Nodes.Dof;
%number of internal nodes: we know that the N unknown nodes are numbered from
%1 to N in Me.UnknownNodes; the maximum is therefore the number of unknown
%(degrees of freedom)
numDof = max(Dof);
spmd(NPools)
    eTr=codistributed(1:size(V, 1));
    
%     HowMany=ceil(size(V, 1)/numlabs);
%     emin=(labindex-1)*HowMany+1;
%     if labindex==numlabs
%         emax=size(V, 1);
%     else
%         emax=labindex*HowMany;
%     end
    %vectors preallocation: instead of allocating the (sparse) diffusion matrix, 
    %we save the rows, columns and values corresponding to each contribution; 
    %at the end, we'll call sparse(...) to obtain the diffusion matrix
    b = zeros(numDof, 1);
    row = zeros(Me.MatrixContributions, 1);
    col = zeros(Me.MatrixContributions, 1);
    d = zeros(Me.MatrixContributions, 1);
    pos = 1;  %we start from the element in position 1, we'll increase this index 
            %everytime we add an entry

    %main loop on each triangle        
    for e = drange(1:size(V,1)),
        Dx(1) = Nodes.X(V(e, 3)) - Nodes.X(V(e, 2));
        Dx(2) = Nodes.X(V(e, 1)) - Nodes.X(V(e, 3));
        Dx(3) = Nodes.X(V(e, 2)) - Nodes.X(V(e, 1));
        Dy(1) = Nodes.Y(V(e, 3)) - Nodes.Y(V(e, 2));
        Dy(2) = Nodes.Y(V(e, 1)) - Nodes.Y(V(e, 3));
        Dy(3) = Nodes.Y(V(e, 2)) - Nodes.Y(V(e, 1));

        Area=Areas(e);
        %we evaluate the external force in the center of mass of this triangle
        force = f(CenterOfMass.X(e),CenterOfMass.Y(e));
        %evaluate the value of the coefficient in front of the Laplace operator
        c = Me.evaluateProperty('c', e);
        %for each vertex of this triangle 
        for ni = 1:3
            %look at the "unknown" numbering: if the node is positive, it
            %corresponds to a degree of freedom of the problem
            ii = Dof(V(e, ni));
            %is it unknown?
            if ii > 0 
                dtmp=c * (Dy(ni) * Dy(ni) + Dx(ni) * Dx(ni)) / (4 * Area);
                row(pos) =ii;
                col(pos) = ii;
                d(pos) = dtmp;
                pos=pos+1;
                %yes it is! second loop on the vertices
                for nj = 1:3
                    jj = Dof(V(e, nj));
                    %%is it unknown as well?
                    if jj > ii
                        %add the contribution to the stiffness matrix 
                        dtmp=c * (Dy(ni) * Dy(nj) + Dx(ni) * Dx(nj)) / (4 * Area);
                        row(pos) =ii;
                        col(pos) = jj;
                        d(pos) = dtmp;
                        pos = pos + 1;
                        row(pos) =jj;
                        col(pos) = ii;
                        d(pos) = dtmp;
                        pos = pos + 1;
                        %Non sparse solution: D(ii,jj)=D(ii,jj) + c*(Dy(i)*Dy(j)+Dx(i)*Dx(j))/(4.0*Area) ;
                    end
                end
                %build the constant terms vector adding the external
                %contribution
                b(ii) = b(ii) + Area * force / 3;
            end
        end
    end
    %assemble the stiffness matrix D
    D = sparse(row(1:pos-1), col(1:pos-1), d(1:pos-1), numDof, numDof);
end
Dfinal=spalloc(numDof,numDof,Me.MatrixContributions);
bfinal=zeros(numDof,1);
for k=1:numel(D)
    Dfinal=Dfinal+D{k};
    bfinal=bfinal+b{k};
end
