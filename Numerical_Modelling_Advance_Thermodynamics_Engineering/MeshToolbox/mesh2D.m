classdef mesh2D
%Class to create and handle 2D meshes

%Version 2017.1
%Copyright 2014-2017 Paolo Bardella
properties (SetAccess=protected)
    %3xNt: one row for each triangle; on each row the index (global numbering) 
    %of the vertices of the triangle
    Triangles               
    %2xNv: one row for each vertex; 
    %first column: x-coordinates, second column: y-coordinates
    Nodes %Coordinates             
    %1xNv: for each vertex, a value >0 indicates that the vertex is a 
    %degree of freedom (and the value is the number of the node in the 
    %unknown vertices numbering); otherwise, the modulus of the (negative)
    %value if the index of row in DirichletNodes for that vertex
    Edges
    BC    
    MatrixContributions
    Time
    Regions
end
properties (Access=protected)
    TriangleFields
end
methods
    function obj = mesh2D(varargin)
    %mesh2D creates and returns a mesh2d object.
    %Me=mesh2D(Region) creates a 2d mesh of the supplied region(s) object(s) using
    %   a maximum traingle size whih is 1/10 of the shortest side of the
    %   rectangle enclosing the region(s)
    %Me=mesh2D(Region, Refine) uses the second parameter as maximum size of
    %   the triangle. Refine can be a scalar, or a function of two
    %   parameters (x and y): in the latter case the size is estimated in
    %   each node of the domain, to allow for selective refinements.
    %Me=mesh2D(Region, Refine, Smooting) applies an additional smoothing to
    %the mesh
    %Me=mesh2D(Region, OldMesh) can be used to efficiently generate a mesh
    %starting from another meshalready calculated on the same region (but
    %with different B.C.s). In this way, since the domain is the same, the
    %nodes are placed in the same position, and only the new B.C.s are
    %taken into account.
        if nargin==0
            error('mesh2D:mesh2D:InvalidInputParameters',...
                'This function requires at least a Region object');
        end
        MaximumArea=[];
        RefinementFunction=[];
        if isa(varargin{1},'region')
            Sh=varargin{1};
        else
            error('mesh2D:mesh2D:InvalidInputParameters',...
                'The first parameter must be a Region object');
        end
        
        if nargin>1
            tmp=varargin{2};
            if isscalar(tmp) && isreal(tmp) && tmp>0
                MaximumArea=tmp;
            elseif isa(tmp,'function_handle') && nargin(tmp)==2
                RefinementFunction=tmp;
            end
        end
        if nargin>2
            error('mesh2D:mesh2D:InvalidInputParameters',...
                'Too many input arguments');
        end
            obj=mesh2D.CreateTriangleMesh(obj,Sh, MaximumArea,RefinementFunction);
%         elseif isa(varargin{1},'mesh2D')
%             MeOld=varargin{1};
%             if nargin==1 || ~isa(varargin{2},'region')
%                 error('mesh2D:mesh2D:InvalidInputParameters',...
%                     'When calling mesh2D passing a mesh2D object as first parameter, the second parameter must be either a scalar or a Region object');
%             end
%             Sh=varargin{2};
%             if nargin>2
%                 error('mesh2D:mesh2D:InvalidInputParameters',...
%                     'Too many input arguments');
%             end
%             obj=mesh2D.UpdateTriangleMesh(obj, MeOld,Sh, MaximumArea);
        
    end        
    
    function uu=copyToAllNodes(Me,u,uu)
    %CopyToAllNodes copies the values u calculated in the unknown nodes
    %in a vector whose size is equal to the total number of nodes of the Mesh
    %The resulting vector can be passed to function Draw() to plot the result
    %If the third parameter is provided, it must be a column vector  with 
    %a number of rows equal to the number of mesh nodes. In this case the
    %values are copied directly in uu.
        MyNarginchk(nargin,2,3);
        if nargin==2
            uu=zeros(size(Me.Nodes.Dof));
        else            
            if any(size(Me.Nodes.Dof)~=size(uu))
            error('mesh2D:mesh2D:InvalidInputParameters',...
            'the third parameter must be a column vector with a number of rows equal to the number of mesh nodes'); 
            end
        end

        uu(Me.Nodes.Dof>0)=u;
        if ~isempty(Me.BC.DirichletNodes)
            uu(Me.Nodes.Dof<0)=Me.BC.DirichletNodes(:,2);
        end
    end
    
    function res=evaluateProperty (Me,propertyName,e)
        %a bit slow, I removed the check for better performance
        %MyNarginchk(nargin,2,3);         
        com=Me.Triangles.CenterOfMass;
        if nargin==2 
            e=1:numel(Me.Triangles.Region);
        end
        
        if isscalar(e)
            res=Me.Regions(Me.Triangles.Region(e)).evaluateProperty(propertyName,com.X(e), com.Y(e));    
        else
        
        %identify the different regions corresponding to the passed
        %triangles indices
        RegionsIndices=unique(Me.Triangles.Region(e)).';
        Prop{length(RegionsIndices)}=[];
        IsNumProp=false(size(RegionsIndices));
        for ks=RegionsIndices
            Prop{ks}=Me.Regions(ks).Properties.(propertyName);
            if isnumeric(Prop{ks})
                IsNumProp(ks)=true;
                if ks==1
                    DataSize=size(Prop{ks},2);
                end
            else
                if ks==1
                    DataSize=size(Prop{ks}(1,1),2);
                end
            end                
        end
        res=zeros(numel(e),DataSize);

        for ks=RegionsIndices
            ThisRegionTriangles=Me.Triangles.Region(e)==ks;
            if IsNumProp(ks)
                if isscalar(Prop{ks})
                    res(ThisRegionTriangles)=Prop{ks};
                else
                    NumTriangles=sum(ThisRegionTriangles);
                    if size(Prop{ks},1)==NumTriangles
                        res(ThisRegionTriangles,:)=Prop{ks};
                    elseif size(Prop{ks},1)==1
                        res(ThisRegionTriangles,:)=repmat(Prop{ks},NumTriangles,1);
                    else
                        error('mesh2D:evaluateProperty:incompatiblesizes',...
                        'The supplied parameter size is not compatible with the size of the triangles');                    
                    end
                end
            else
                res(ThisRegionTriangles,:)=Prop{ks}(com.X(ThisRegionTriangles), ...
                    com.Y(ThisRegionTriangles));
            end
        end
        end
    end
    
    function [dx,dy]=gradient(Me,uu)
    %gradient returns the gradients in the x- and y- directions of the
    %solution uu (calculated in ALL the nodes)    
        MyNarginchk(nargin,2,2);
        V=Me.Triangles.Vertices;
        A=Me.Triangles.Areas;
        N=Me.Nodes;
        dx=zeros(size(Me.Nodes.X,1),1);
        dy=zeros(size(Me.Nodes.X,1),1);
        freq=zeros(size(Me.Nodes.X,1),1);
        %I evaluate the gradients on each triangle
        for e=1:size(V,1)
            Dx(1) = N.X(V(e,3)) - N.X(V(e,2));
            Dx(2) = N.X(V(e,1)) - N.X(V(e,3));
            Dx(3) = N.X(V(e,2)) - N.X(V(e,1));

            Dy(1) = -N.Y(V(e,3)) + N.Y(V(e,2));
            Dy(2) = -N.Y(V(e,1)) + N.Y(V(e,3));
            Dy(3) = -N.Y(V(e,2)) + N.Y(V(e,1));

            Dz(1) = uu(V(e,2)) - uu(V(e,3));
            Dz(2) = uu(V(e,3)) - uu(V(e,1));
            Dz(3) = uu(V(e,1)) - uu(V(e,2));

            gry=-1/(2*A(e))*(Dx(3)*Dz(2)-Dx(2)*Dz(3));
            grx=1/(2*A(e))*(Dy(3)*Dz(2)-Dz(3)*Dy(2));
            dx(V(e,:))=dx(V(e,:))+[grx;grx;grx];
            dy(V(e,:))=dy(V(e,:))+[gry;gry;gry];
            freq(V(e,:))=freq(V(e,:))+1;
        end
        dx=dx./freq;
        dy=-dy./freq;
    end
    function I=integrate(Me,uu,Region)
    %integrate returns the integral solution uu 
    %I=integrate(Me,uu) returns a scalar indicating the numerical integral
    %calculated over the whole domain
    %I=integrate(Me,uu, Region) returns the integral calculated over each 
    %specified Region
    
        MyNarginchk(nargin,2,3);
        if nargin==2
            Region=unique(Me.Triangles.Region);
        end
        Tr=Me.Triangles.Vertices;
        A=Me.Triangles.Areas;
        I=zeros(numel(Region),1);
        for f=1:numel(Region)
            for e=1:length(Tr)
                if(Me.Triangles.Region(e)==Region(f))
                    I(f)=I(f)+mean(uu(Tr(e,:)))*A(e);
                end
            end
        end
    end
    
    function z=interpolate(Me, uu, xy, tr)
        if nargin==3
            [~, tr]=findClosestNode(Me,xy);
        end
        z=zeros(size(xy,1),1);
        V=Me.Triangles.Vertices;
        C=Me.Nodes;
        for kz=1:length(z)
            e=tr(kz);
            if e==-1
                z(kz)=NaN;
            else
                x= C.X(V(e,:));  
                y= C.Y(V(e,:));    
                u=uu(V(e,:));
                dx2 = x(2) - x(1);    dy2 = y(2) - y(1);
                dx3 = x(3) - x(1);    dy3 = y(3) - y(1);    
                dx0 = xy(kz,1)-x(1);  dy0 = xy(kz,2)-y(1);
                dz2 = u(2) - u(1);    dz3 = u(3) - u(1);    
                z(kz)=u(1)+(-dx0*(dy2*dz3-dy3*dz2)+dy0*(dx2*dz3-dx3*dz2))/(dx2*dy3-dx3*dy2);
            end
        end
    end
    function [PosNodes,TriangleIndices]=findClosestNode(Me,xyv)
    %findClosestNode returns the index, in the global numbering, of the
    %node which is closest to the assigned xyv node.
    %Input:
    %   Me:     mesh2d object
    %   xyv:    vector or matrix with 2 columns (x and y coordinates)
    %Output:
    %   PosNodes:   indices of the closest nodes to the point specified in 
    %               the corresponding row of xyv.
    %   TriangleIndices: indices of the triangles in which the points in
    %               xyv are included. The value is -1 if the point is
    %               outside the domain
    %Example:
    %Me=mesh2D(regions.rect());
    %Me.findClosestNode([0,0.1; 0.2,.3]);
        if nargin~=2
            error('mesh2D:findClosestNode:InvalidInputParameters',...
                'The function requires 2 parameters');
        end
        if size(xyv,2)~=2 || ~isreal(xyv) ||~isnumeric(xyv)
            error('mesh2D:findClosestNode:InvalidInputParameters',...
                'The second parameter should be a 2 columns real matrix');
        end
        PosNodes=zeros(size(xyv,1),1);       
        X=Me.Nodes.X;
        Y=Me.Nodes.Y;
        for kxyv=1:size(xyv,1)
            xy=xyv(kxyv,:);
            [~, PosNodes(kxyv)]=min((X-xy(1)).^2+(Y-xy(2)).^2);
        end
        if nargout==2
            TriangleIndices=zeros(size(xyv,1),1);   
            %calculate circumcenters
            S2 = X.^2+Y.^2;
            V=Me.Triangles.Vertices;
            Areas=Me.Triangles.Areas;
            X1=X(V(:,1)); X2=X(V(:,2)); X3=X(V(:,3));
            Y1=Y(V(:,1)); Y2=Y(V(:,2)); Y3=Y(V(:,3));
            d12X=X1-X2;    d12Y=Y1-Y2;
            d23X=X2-X3;    d23Y=Y2-Y3;
            d31X=X3-X1;    d31Y=Y3-Y1;
            D  = 0.5 ./ (X1.*d23Y + X2.*d31Y + X3.*d12Y);
            % center coordinates
            xc =  (S2(V(:,1)).*d23Y + S2(V(:,2)).*d31Y + S2(V(:,3)).*d12Y ).*D;
            yc = -(S2(V(:,1)).*d23X + S2(V(:,2)).*d31X + S2(V(:,3)).*d12X ).*D;
            r2=(xc-X1).^2+(yc-Y1).^2; %radius squared
            for kxyv=1:size(xyv,1)
                xn=xyv(kxyv,1);
                yn=xyv(kxyv,2);
                %for the selected point, find the triangles for which the
                %distance p-(xc,yc) is less than the circumradius
                guessTr=find((xn-xc).^2+(yc-yc).^2<=r2);
                %guessTr are the indice of the possible "good" triangles
                %for each one check if the node falls inside
                if isempty(guessTr)
                    TriangleIndices(kxyv)=-1;
                else
                    err=zeros(size(guessTr));
                    for k=1:numel(guessTr)
                        x1=X1(guessTr(k));      y1=Y1(guessTr(k));
                        x2=X2(guessTr(k));      y2=Y2(guessTr(k));
                        x3=X3(guessTr(k));      y3=Y3(guessTr(k));
                        err(k)=Areas(guessTr(k))-0.5*...
                            (abs(xn*(y1-y2)-x1*(yn-y2)+x2*(yn-y1))+...
                             abs(xn*(y1-y3)-x1*(yn-y3)+x3*(yn-y1))+...
                             abs(xn*(y2-y3)-x2*(yn-y3)+x3*(yn-y2)));
                        
                           % (0.5*abs(det([xn, x1, x2; yn, y1, y2;1,1,1]))+... %A(xy,xy1,xy2)-...
                           %  0.5*abs(det([xn, x1, x3; yn, y1, y3;1,1,1]))+... %A(xy,xy1,xy3)-...
                           %  0.5*abs(det([xn, x2, x3; yn, y2, y3;1,1,1]))); %A(xy,xy2,xy3)-...
                    end
                    [val,pos]=min(abs(err)./eps(Areas(guessTr)));
                     %troppo grande, non appartiene ai triangoli individuati
                    if val>120
                        TriangleIndices(kxyv)=-1;
                    else
                        TriangleIndices(kxyv)=guessTr(pos);
                    end
                end
            end                
        end
    end
        
    function Result=find(Me,Property, Kind)
    %find returns a boolean vector indicating, for each node of the mesh if
    %function Property is true in that node.
    %Result=find(Me,Property) where Me is a mesh2d object and Property is a
    %function handle accepting 2 inputs (x and y coordinates) evaluates Property 
    %for each node in the grid and returns the boolean result of the evaluation
    %Result=find(Me,Property, Kind) , where kind can be either 
    %'dirichlet','unknwon', or 'all', evaluate Property only to the
    %Dirichlet, or Unknown nodes, or to all the nodes (the latter is the
    %default behaviour).
    %
    %Example
    %Me=mesh2D(regions.rect());
    %Me.find(@(x,y)x<y,'dirichlet') %returns the indices of the dirichlet 
    %nodes having x-coordinate < y-coordinate
    
        if nargin<2
            error('mesh2D:find:WrongNumberInputParameters',...
                'The function requires 2 at least parameters');
        end
        if nargin<3
            Kind='all';
        end
        if ~isa(Property,'function_handle') 
            error('mesh2D:find:InvalidInputParameters',...
                'The second parameter should be a function with 2 input arguments');
        end
        if nargin(Property) ~=2
            error('mesh2D:find:InvalidInputParameters',...
                'The second parameter should be a function with 2 input arguments');
        end
        try
            Result=Property(Me.Nodes.X,Me.Nodes.Y);
        catch Exception
            disp(Exception);
            error('mesh2D:find:WrongProperties',...
                'The evaluation of the supplied function returned an error');
        end
        if strcmpi(Kind,'dirichlet')==1 ||strcmpi (Kind,'d')==1
            Result=(Result & Me.Nodes.Dof<0);
        elseif strcmpi(Kind,'unknown')==1 ||strcmpi (Kind,'u')==1
            Result=(Result & Me.Nodes.Dof>0);
        elseif strcmpi(Kind,'all')==0 && strcmpi (Kind,'a')==0
            error('mesh2D:find:WrongKind',...
                'The kind of the returned nodes must be either dirichlet, unknown, or all');
        else
            Result=(Result & 1); %force conversion to boolean
        end            
    end
    
    function obj=forceDirichlet(obj, Node, DirichletValue)
    %forceDirichlet convert an unknown node to a Dirichlet condition.
    %forceDirichlet(Me, Node, DirichletValue) assigns a Dirichlet condition 
    %to the node Node in the mesh object. The value of the Dirichlet 
    %condition in the node is inditted by the third parameter, DirichletValue
    
        MyNarginchk(nargin,3,3);
        if Node<=0 || Node>size(obj.Nodes.Dof,1)
           error('mesh2D:forceDirichlet:WrongInput',...
               ['The Node must be between 1 and ' num2str(size(obj.Nodes.Dof,1))]);
        end
        if obj.Nodes.Dof(Node)<0 %already dirichlet->update the value
            obj.BC.DirichletNodes(-obj.Nodes.Dof(Node),2)=DirichletValue;
        else
            obj.BC.DirichletNodes=[obj.BC.DirichletNodes;[Node,DirichletValue]];
            obj.Nodes.Dof(Node)=-size(obj.BC.DirichletNodes,1);
            obj.Nodes.Dof(obj.Nodes.Dof>0)=1:sum(obj.Nodes.Dof>0);
            obj.MatrixContributions=sum(sum(obj.Nodes.Dof(obj.Triangles.Vertices)>0,2).^2);
        end
        if nargout==0
            warning('mesh2D:forceDirichlet:ResultNotSaved', ...
                'The node has been set, but the result was not saved in a variable');
        end
    end
    
    function objNew=extractMesh(obj,RegionsToKeep)
        RegionsToKeep=unique(RegionsToKeep);
        if RegionsToKeep(end)>length(obj.Regions)
            error('You asked for region %d but your mesh has only %d region(s)',RegionsToKeep(end), length(obj.Regions));
        end
        TotalNewRegions=length(RegionsToKeep);
        TotalOriginalRegions=length(obj.Regions);
        Region_renumbering=zeros(TotalOriginalRegions,1);
        Region_renumbering(RegionsToKeep)=1:TotalNewRegions;
        objNew=obj;
        TotalOriginalNodes=length(objNew.Nodes.X);        
        %first update the Triangles. Only keep the ones corresponging to
        %the selected region
        TrianglesToKeep=ismember(obj.Triangles.Region,RegionsToKeep);
        TotalNewTriangles=nnz(TrianglesToKeep);
        objNew.Triangles.Areas=obj.Triangles.Areas(TrianglesToKeep);
        objNew.Triangles.CenterOfMass.X=obj.Triangles.CenterOfMass.X(TrianglesToKeep);
        objNew.Triangles.CenterOfMass.Y=obj.Triangles.CenterOfMass.Y(TrianglesToKeep);
        objNew.Triangles.Vertices=obj.Triangles.Vertices(TrianglesToKeep,:);
        objNew.Triangles.Region=Region_renumbering(objNew.Triangles.Region(TrianglesToKeep));
        %Which nodes should we keep?
        NodesToKeep=unique(objNew.Triangles.Vertices(:));
        TotalNewNodes=length(NodesToKeep);
        Node_renumbering=zeros(TotalOriginalNodes,1);
        Node_renumbering(NodesToKeep)=1:TotalNewNodes;
        objNew.Triangles.Vertices=Node_renumbering(objNew.Triangles.Vertices);
        objNew.Triangles.Vertices=reshape(objNew.Triangles.Vertices,TotalNewTriangles,3);
        
        
        %update the nodes property
        objNew.Nodes.X=obj.Nodes.X(NodesToKeep);
        objNew.Nodes.Y=obj.Nodes.Y(NodesToKeep);
        objNew.Nodes.Dof=obj.Nodes.Dof(NodesToKeep);
        
        %which edges should we keep? only the ones connecting 2 kept nodes.
        EdgesToKeep=all(reshape(ismember(objNew.Edges(:),NodesToKeep),[],2).').';
        IndicesEdgesToKeep=find(EdgesToKeep);
        TotalNewEdges=nnz(EdgesToKeep);
        
        TotalOriginalEdges=size(obj.Edges,1);        
        Edge_renumbering=zeros(TotalOriginalEdges,1);
        Edge_renumbering(EdgesToKeep)=1:TotalNewEdges;        
        
        objNew.Edges=objNew.Edges(EdgesToKeep,:);
        objNew.Edges=Node_renumbering(objNew.Edges);
        objNew.Edges=reshape(objNew.Edges,TotalNewEdges,2);
        
        NumDofNew=length(objNew.Nodes.Dof(objNew.Nodes.Dof>0));
        NNonDofNew=TotalNewNodes-NumDofNew;

        objNew.BC.DirichletNodes=obj.BC.DirichletNodes(-objNew.Nodes.Dof(objNew.Nodes.Dof<0),:);

        objNew.BC.DirichletNodes(:,1)=Node_renumbering(objNew.BC.DirichletNodes(:,1));
        objNew.BC.NeumannEdges=objNew.BC.NeumannEdges(ismember(objNew.BC.NeumannEdges(:,1),IndicesEdgesToKeep) ,:);
        objNew.BC.NeumannEdges(:,1)=Edge_renumbering(objNew.BC.NeumannEdges(:,1));
        
        objNew.BC.RobinEdges=objNew.BC.RobinEdges(ismember(objNew.BC.RobinEdges(:,1),IndicesEdgesToKeep) ,:);
        objNew.BC.RobinEdges(:,1)=Edge_renumbering(objNew.BC.RobinEdges(:,1));
        
        objNew.Nodes.TwinNode(objNew.Nodes.TwinNode>0)=Node_renumbering(objNew.Nodes.TwinNode(objNew.Nodes.TwinNode>0));
        objNew.Nodes.TwinNode=objNew.Nodes.TwinNode(NodesToKeep);       

        objNew.Nodes.Dof(objNew.Nodes.Dof>0)=1:NumDofNew;
        objNew.Nodes.Dof(objNew.Nodes.Dof<0)=-(1:NNonDofNew);
        objNew.Regions=objNew.Regions(RegionsToKeep);
        objNew.MatrixContributions=sum(sum(objNew.Nodes.Dof(objNew.Triangles.Vertices)>0,2).^2);
        
        %Update NeumannEdges to include all the edges that were previously
        %in contact with other regions
        SortedEdges=sort(objNew.Edges,2);
        TrSide=zeros(size(objNew.Triangles.Vertices));
        for k=1:size(objNew.Triangles.Vertices,1)
            v=sort(objNew.Triangles.Vertices(k,:));
            TrSide(k,1)=find(SortedEdges(:,1)==v(1) & SortedEdges(:,2)==v(2));
            TrSide(k,2)=find(SortedEdges(:,1)==v(2) & SortedEdges(:,2)==v(3));
            TrSide(k,3)=find(SortedEdges(:,1)==v(1) & SortedEdges(:,2)==v(3));
        end
        BorderEdges=find(hist(TrSide(:),unique(TrSide(:)))==1);
        NeumannOrRobinBorders=BorderEdges((any(not(ismember(SortedEdges(BorderEdges,:),objNew.BC.DirichletNodes(:,1))'))));
        %Now, remove all the edges which are already listed as Neumann or
        %Robin
        AlreadyDefinedBorders=[objNew.BC.RobinEdges(:,1); objNew.BC.NeumannEdges(:,1)];
        NewNeumann=BorderEdges((not(ismember(NeumannOrRobinBorders,AlreadyDefinedBorders)))).';
        objNew.BC.NeumannEdges=[objNew.BC.NeumannEdges;[NewNeumann,zeros(size(NewNeumann))]];
    end

    function draw(varargin)
        %DrawMesh(Me) receives a structure from createMesh and draws it
        %Possible additional options are
        % e or edge or edges:
        % n or node or nodes:
        % t or triangle or triangles:
        % i or internalnode or internalnodes:
        % d or dirichlet or dirichletnodes or dirichletnode:
        % q or trianglequality or quality
        % c or tricontour or contour
        % o or offset
        % h or hidemesh
        % s or showmesh 
        % p or periodic
        Me=varargin{1};
        X=Me.Nodes.X;
        Y=Me.Nodes.Y;
        t=Me.Triangles.Vertices;
        dof=Me.Nodes.Dof;
        com=Me.Triangles.CenterOfMass;
        offset=zeros(size(X));
        fnum=Me.Triangles.Region;
        Opt=mesh2D.drawMeshOptions(varargin(2:end));
        if ~isempty(Opt.Solution)
            if Opt.TriContour
                tricontour([X, Y],t,Opt.Solution,20);
            else
                trisurf(t,X,Y,Opt.Solution);
                if Opt.Offset
                    offset=offset+min(Opt.Solution);
                end
            end
            colorbar;
            hold on;
        end
        face=unique(fnum);
        
% % % %         if Opt.TriangleQuality && ~isempty(Me.Triangles.Quality)
% % % %             quality=Me.Triangles.Quality;
% % % %         elseif kend<length(varargin)
% % % %             name=varargin{kend};
% % % %             quality=zeros(length(Me.Triangles),1);
% % % %             for e=1:length(Me.Triangles)
% % % %                 xG = Me.CenterOfMass(e,1);
% % % %                 yG = Me.CenterOfMass(e,2);
% % % %                 s=fnum(e);
% % % %                 quality(e)= get(Me.Regions(s),name,xG,yG);
% % % %             end
% % % %         end 
%         % Colour mesh for each face
        if Opt.ShowMesh~=0
            %plot3(X,Y,offset,'b.','markersize',1)
            hold on;
            col = ['b','r','g','y','m'];
            for k = 1:length(face)
               
                    colk = mod(k,length(col));
                    if (colk==0)
                        colk = length(col);
                    end
                    patch('faces',t(fnum==face(k),:),'vertices',[X, Y, offset],'facecolor','none','edgecolor',col(colk));
                
            end
        end        
        hold on;
        %% Dirichlet nodes
        if Opt.TextOnDirichletNodes
              ind=Me.BC.DirichletNodes(:,1);
              %%plot Dirichlet nodes on that face              
              plot3(X(ind),Y(ind),offset(ind),['o' col(colk)]);    
        end

        %% Number each node
        if Opt.TextOnNodes
            Num=length(X);
            for n=1:Num,
                text(X(n),Y(n),offset(n),num2str(n) );
            end 
        end

        %% number each triangle
        if Opt.TextOnTriangles
            for kt=1:length(t),
                text(com.X(kt),com.Y(kt),offset(1),num2str(kt),...
                    'horizontalalignment','center','color',[1 0 0]');
                %per centrare orizzontalmente il testo 
            end 
        end

        %% number each edge
        if Opt.TextOnEdges
            for e=1:length(Me.Edges),
                Vert=Me.Edges(e,[1 2]); 		%3 vertices 
                x=mean(X(Vert));	%x average
                y=mean(Y(Vert));	%y average
                %red text, horizontally centered
                text(x,y,offset(1),num2str(e),...
                    'horizontalalignment','center','color',[1 0 0]');                 
            end 
        end

        %% number each TextOnInternalNodes
        if Opt.TextOnInternalNodes
            for n=1:length(dof)
                if(dof(n)>0)
                    text(X(n),Y(n),offset(n),num2str(dof(n)));
                end
            end
        end
        
        %% periodic nodes test
        if Opt.Periodic
            Periodic=find(Me.Nodes.TwinNode>0);            
            plot(Me.Nodes.X(Periodic),Me.Nodes.Y(Periodic),'og');
        end
        %% set view
        if isempty(Opt.Solution)
            view(0,90);
        end
        xlabel('x dir [m]');
        ylabel('y dir [m]');        
    end
end

methods (Access=protected, Static)
    function [in, GlobalConditions]= GenerateInputToTriangle(Sh,MaximumArea)
        in.verboseinput=0;
        if not(isempty(MaximumArea))
            in.maxArea=MaximumArea;
        else
            M=max(Sh);
            m=min(Sh);
            in.maxArea=abs(prod(M-m))/25;
        end
        %holes: sum all the regions
        ShTot=Sh(1);
        for k=2:numel(Sh)
            ShTot=ShTot+Sh(k);
        end
        
        in.numberofholes=sum([ShTot.Borders(:).Hole]==1);        
        in.holelist=zeros(in.numberofholes,2);
        numberofholes=1;
        for kb=1:length(ShTot.Borders)
            if ShTot.Borders(kb).Hole==1
                %add to the holes
                [X,Y]=mesh2D.GetInternalPoint(ShTot.Borders(kb).X,ShTot.Borders(kb).Y,true);        
                in.holelist(numberofholes,:)=[X,Y];
                numberofholes=numberofholes+1;
            end
        end
        
        %integrate the number of nodes on edges with Periodic B.C.s
        %I suppose that triangle does not return angles>90'. The minimum
        %angle is 20'. When we have a 90' angle, the area is 
        %Area=l^2*sin(alpha)/2, but alpha is >20', so
        %lmax=sqrt(2*Area/sin(20'))=2.5*sqrt(A). I then divide by 2 to be
        %conservative
        
        lMax=2.5*sqrt(in.maxArea)/2;
        for ks=1:length(Sh)
            %for each border
           
            for kb=1:numel(Sh(ks).Borders)
                 XPeriodic={};
                YPeriodic={};
                After={};
                TotInserted=0;
                for kbc=1:numel(Sh(ks).Borders(kb).Bc)
                    if Sh(ks).Borders(kb).Bc(kbc).Condition=='P'
                        %calculate the length of this side
                        kStart=kbc;
                        kEnd=rem(kbc,numel(Sh(ks).Borders(kb).Bc))+1;
                        X=Sh(ks).Borders(kb).X;
                        Y=Sh(ks).Borders(kb).Y;
                        dist=sqrt((X(kStart)-X(kEnd))^2+(Y(kStart)-Y(kEnd))^2);
                        %if the side is too long, split it adding nodes 
                        if dist>lMax
                            XPeriodic{end+1}=linspace(X(kStart),X(kEnd),floor(dist/lMax)).';%#ok
                            XPeriodic{end}=XPeriodic{end}(2:end-1);
                            YPeriodic{end+1}=linspace(Y(kStart),Y(kEnd),floor(dist/lMax)).'; %#ok 
                            YPeriodic{end}=YPeriodic{end}(2:end-1);
                            After{end+1}=kStart+TotInserted;%#ok
                            TotInserted=TotInserted+numel(YPeriodic{end});
                        end
                    end
                end
                if TotInserted>0
                    for kA=1:length(After)
                        Sh(ks).Borders(kb).insertNode(After{kA},[XPeriodic{kA},YPeriodic{kA}]);
                    end
                end
            end
        end
        
        %calculate number of nodes and regions
        in.numberofpointattributes=0;
        in.numberofpoints=0;
        in.numberofregions=0;
        %for each region
        for ks=1:length(Sh)
            %for each border
            for kb=1:numel(Sh(ks).Borders)
                %number of nodes
                in.numberofpoints=in.numberofpoints+size(Sh(ks).Borders(kb).X,1);        
                if Sh(ks).Borders(kb).Hole==0 
                    %not an hole
                    in.numberofregions=in.numberofregions+1;   
                end
            end
        end

        %allocate according to the calculated quantities
        in.regionlist= zeros(in.numberofregions,4);
        Coordinates=zeros(in.numberofpoints,2);
        Segments=zeros(in.numberofpoints,2);

        pos=1;
        numberofregions=1;

        for ks=1:length(Sh)    
            if isfield(Sh(ks).Properties,'MeshMaxArea')
                area=Sh(ks).Properties.MeshMaxArea;
                if ~isnumeric(area) || ~isscalar(area) || area<=0
                    error('mesh2D:forceDirichlet:WrongParamater',...
                        'Property MeshMaxArea must be a positive scalar number');
                end
            else
                area=in.maxArea;
            end
                
            for kb=1:numel(Sh(ks).Borders)    
                %fill coordinates and segments
                %number of nodes in this border
                NumLocalNodes=length(Sh(ks).Borders(kb).X);
                %array of positions in coordinates and segments
                Indices=pos:pos+NumLocalNodes-1;
                Coordinates(Indices,:)=[Sh(ks).Borders(kb).X,Sh(ks).Borders(kb).Y];
                %generate two vectors: one goes from form 1 to N, the other one
                %starts with N, and them from 1 to N-1
                Segments(Indices,:)=[Indices',[Indices(end), Indices(1:end-1)]'];
                %save all the BC in an array
                GlobalConditions(Indices)=[Sh(ks).Borders(kb).Bc(end);Sh(ks).Borders(kb).Bc(1:end-1)];%#ok
                %increase the position for the next border/region
                pos=pos+NumLocalNodes;
                %obtain the coordinates of a point inside the region defined by the 
                %current border
                if Sh(ks).Borders(kb).Hole==0
                    [X,Y]=mesh2D.GetInternalPoint(Sh(ks).Borders(kb).X,Sh(ks).Borders(kb).Y,false);        
                    %add to the region.                     
                    in.regionlist(numberofregions,:)=[X;Y;ks;area];            
                    numberofregions=numberofregions+1;
                end
            end
        end


        %remove duplicated coordinates
        [Coordinates,~,ib]=unique(Coordinates,'rows'); 
        %and reshape Segments accordingly
        Segments=reshape(ib(Segments(:)),[],2);
        %remove duplicated Segments. Sort is used to consider segments [a b]
        %equivalent to [b a]
        [Segments,ia,ib]=unique(sort(Segments,2),'rows');

        in.numberofpoints=size(Coordinates,1);
        in.pointlist=Coordinates;
        in.numberofsegments=size(Segments,1);
        in.segmentslist=Segments;
        in.segmentmarkerlist=(1:size(Segments,1))';

        %loop over all the edges
        for k=1:numel(ia), 

            OverlappedEdges=(ib==k); 
            if nnz(OverlappedEdges)>1%we have two or more edges connecting the same 2 nodes
                %should we keep inner Dirichlet BC?
                if exist('DirichletOnSegments','var') && DirichletOnSegments==1
                    if any([GlobalConditions(OverlappedEdges).Condition]=='D')
                        %must check is the conditions are the same.
                        %We cannot check equality for handles or anonymous
                        %functions, so we use this logic: if the BC are all the
                        %same number, then keep Dirichlet, Otherwise, set
                        %Continuity
                        first=find(OverlappedEdges,1);
                         if isempty(GlobalConditions(first).Value) ||...
                            any([GlobalConditions(OverlappedEdges).Value]~=GlobalConditions(first).Value)
                            warning('changing BC');
                            GlobalConditions(OverlappedEdges)=boundaries.continuity;
                         else
                            GlobalConditions(OverlappedEdges)=GlobalConditions(first);
                         end
                    else
                        GlobalConditions(OverlappedEdges)=boundaries.continuity;
                    end
                else
                    %the edge is internal respect to the domain; assign continuity
                    GlobalConditions(OverlappedEdges)=boundaries.continuity;
                end
            end
        end
        %keep the BC for the non duplicated edges only
        GlobalConditions=GlobalConditions(ia);
        %now assign the nodes to an edge/BC
        PointsMarker=zeros(size(Coordinates,1),1);
        %for each node
        for k=1:size(Coordinates,1)
            %find to which edges it belongs to
            SegmentsWithNodek=find(Segments(:,1)==k | Segments(:,2)==k);
            %and the corresponding BCs
            Condition=[GlobalConditions(SegmentsWithNodek).Condition];
            %is there at least one edge with Dirichlet? 
            pos=find(Condition=='D',1);
            if ~isempty(pos)
                %if so, assign the first Dirichlet edge. 
                PointsMarker(k)=SegmentsWithNodek(pos);
            else
                %is there at least one edge with periodic condition? 
                pos=find(Condition=='P',1);
                if ~isempty(pos)
                    PointsMarker(k)=SegmentsWithNodek(pos);
                else
                    %is there at least one edge with Robin condition? 
                    pos=find(Condition=='R',1);
                    if ~isempty(pos)
                        PointsMarker(k)=SegmentsWithNodek(pos);
                    else
                    %is there at least one edge with Robin condition? 
                        pos=find(Condition=='N',1);
                        if ~isempty(pos)
                            PointsMarker(k)=SegmentsWithNodek(pos);
                        else
                            pos=find(Condition=='C',1);
                            if ~isempty(pos)
                                PointsMarker(k)=SegmentsWithNodek(pos);
                            else
                               error('mesh2D:mesh2D:InvalidBoundaryCondition',...
                                    'An edge has an invalid boundary condition');
                            end
                        end
                    end
                end
            end
        end
        in.pointmarkerlist=PointsMarker;        
    end
    function obj=AssignBCs(obj, TFields, GlobalConditions)
        NumNodes=length(obj.Nodes.X);
        GlobalConditionsFrom1=[' ',[GlobalConditions(:).Condition]]; 
        IsDof=GlobalConditionsFrom1(TFields.pointmarkerlist+1)~='D';
        TotDof=sum(IsDof);
        if TotDof==0
            error('mesh2D:mesh2D:InvalidInputParameters',...
                'This geometry has no DOF. Try reducing the maximum triangle area');
        end
        obj.Nodes.Dof=zeros(NumNodes,1);
        obj.Nodes.Dof(IsDof)=1:TotDof;
        obj.Nodes.Dof(~IsDof)=-(1:(NumNodes-TotDof));
        obj.BC.DirichletNodes=[find(~IsDof)', zeros(NumNodes-TotDof,1)];
        DirichletBoundaries=GlobalConditions(TFields.pointmarkerlist(obj.BC.DirichletNodes(:,1)));
        for k=1:length(DirichletBoundaries)
            if isempty(DirichletBoundaries(k).Value)
                xnode=obj.Nodes.X(obj.BC.DirichletNodes(k,1));
                ynode=obj.Nodes.Y(obj.BC.DirichletNodes(k,1));
                obj.BC.DirichletNodes(k,2)=DirichletBoundaries(k).Function(xnode,ynode);
            else
                obj.BC.DirichletNodes(k,2)=DirichletBoundaries(k).Value;
            end
        end
        
        obj.Nodes.TwinNode=zeros(NumNodes,1);
        PeriodicMarkers=find([GlobalConditions(:).Condition]=='P');
        for k=PeriodicMarkers
            Law=GlobalConditions(k);
            for kp=find(TFields.pointmarkerlist==k).',
                TwinXY=Law.Function(obj.Nodes.X(kp),obj.Nodes.Y(kp));                
                [~, posTwin]=min((obj.Nodes.X-TwinXY(:,1)).^2+(obj.Nodes.Y-TwinXY(:,2)).^2);
                obj.Nodes.TwinNode(kp)=posTwin;                
            end
        end
        PeriodicNodes=find(obj.Nodes.TwinNode>0);
        if not(all(obj.Nodes.TwinNode(obj.Nodes.TwinNode(PeriodicNodes))==PeriodicNodes))
            warning('mesh2D:CreateTriangleMesh:InconsistenNumbering', ...
                'One or more nodes with Periodic B.C.s have no twin node.');
        end
        obj.Nodes.TwinNode(PeriodicNodes)=obj.Nodes.Dof(obj.Nodes.TwinNode(PeriodicNodes));

        BorderEdges=find(TFields.edgemarkerlist>0);
        NeumannEdges=BorderEdges([GlobalConditions(TFields.edgemarkerlist(BorderEdges)).Condition]=='N');
        obj.BC.NeumannEdges=[NeumannEdges, zeros(size(NeumannEdges,1),1)];
        NeumannBoundaries=GlobalConditions(TFields.edgemarkerlist(NeumannEdges));
        for k=1:length(NeumannEdges)
            if isempty(NeumannBoundaries(k).Value)
                nodes=obj.Edges(obj.BC.NeumannEdges(k,1),:);                
                xr=mean(obj.Nodes.X(nodes));
                yr=mean(obj.Nodes.Y(nodes));
                h=NeumannBoundaries(k).Function(xr,yr);
                obj.BC.NeumannEdges(k,2)=h;                
            else
                obj.BC.NeumannEdges(k,2)=NeumannBoundaries(k).Value;
            end
        end


        RobinEdges=BorderEdges([GlobalConditions(TFields.edgemarkerlist(BorderEdges)).Condition]=='R');
        obj.BC.RobinEdges=[RobinEdges, zeros(size(RobinEdges,1),2)];
        RobinBoundaries=GlobalConditions(TFields.edgemarkerlist(RobinEdges));
        for k=1:length(RobinEdges)
            if isempty(RobinBoundaries(k).Value)
                nodes=obj.Edges(obj.BC.RobinEdges(k,1),:);                
                xr=mean(obj.Nodes.X(nodes));
                yr=mean(obj.Nodes.Y(nodes));
                h=RobinBoundaries(k).Function{1}(xr,yr);
                g=RobinBoundaries(k).Function{2}(xr,yr);
                obj.BC.RobinEdges(k,[2,3])=[h,g];
            else
                obj.BC.RobinEdges(k,[2,3])=RobinBoundaries(k).Value;
            end
        end


        obj.MatrixContributions=sum(sum(obj.Nodes.Dof(obj.Triangles.Vertices)>0,2).^2);
        if obj.MatrixContributions==0
            error('mesh2D:mesh2D:InvalidInputParameters',...
                'This geometry has no DOF. Try reducing the maximum triangle area');
        end
    end
    function obj=UseTriangleResult(obj, Sh, TriangleStruct)
        obj.Regions=Sh;%original regions  
        
        obj.Nodes.X=TriangleStruct.pointlist(:,1);
        obj.Nodes.Y=TriangleStruct.pointlist(:,2);
        obj.Triangles.Vertices=TriangleStruct.trianglelist;
        x=obj.Nodes.X(obj.Triangles.Vertices);
        y=obj.Nodes.Y(obj.Triangles.Vertices);
        obj.Triangles.CenterOfMass.X=mean(x,2);
        obj.Triangles.CenterOfMass.Y=mean(y,2);
        obj.Triangles.Areas=0.5*abs((x(:,3)-x(:,1)).*(y(:,2)-y(:,1))-(x(:,2)-x(:,1)).*(y(:,3)-y(:,1)));
        obj.Triangles.Region=TriangleStruct.triangleattributelist;
        obj.Edges=TriangleStruct.edgelist;
       
        obj.TriangleFields.pointmarkerlist=sparse(TriangleStruct.pointmarkerlist);    
        obj.TriangleFields.edgemarkerlist=sparse(TriangleStruct.edgemarkerlist);    
    end
    
    function obj=CreateTriangleMesh(obj,Sh, MaximumArea, RefinementFunction)
        [in,GlobalConditions]= mesh2D.GenerateInputToTriangle(Sh, MaximumArea);        
        AbsoluteTime(1)=cputime();        
        if isempty(RefinementFunction)
            TriangleStruct=TriangleMEX(in);
        else
            TriangleStruct=TriangleMEX(in,RefinementFunction);
        end
        %TriangleStruct.GlobalConditions=GlobalConditions;
        AbsoluteTime(2)=cputime();
        obj=mesh2D.UseTriangleResult(obj,Sh, TriangleStruct);
        obj=mesh2D.AssignBCs(obj,obj.TriangleFields,GlobalConditions);
        
        AbsoluteTime(3)=cputime;
        obj.Time=diff(AbsoluteTime);
    end
    
    

%     function obj=UpdateTriangleMesh(obj, objSource,Sh, MaximumArea)
%         [~,GlobalConditions]= mesh2D.GenerateInputToTriangle(Sh, MaximumArea);
%         AbsoluteTime(1)=0;
%         AbsoluteTime(2)=cputime();
%         obj.Triangles=objSource.Triangles;   
%         obj.Nodes=objSource.Nodes;
%         obj.Edges=objSource.Edges;
%         obj.Regions=objSource.Regions;
%         obj=mesh2D.AssignBCs(obj, objSource.TriangleFields,GlobalConditions);
%         AbsoluteTime(3)=cputime;
%         obj.Time=diff(AbsoluteTime);
%     end

    function DMO=drawMeshOptions(OptionsList)
        %Create a struct with the drawing parameters for Draw
        %   Available parameters are:
        %   Edge:       print the number of each edge
        %   Node:       print the number of each node
        %   Triangle:   print the number of each triangle
        %   InternalNode:   print the number of each internal (unknown) node
        %   Dirichletnode:  print the number of each Dirichlet (known) node
        %   Periodic:    highligths the periodic nodes
        DMO.TextOnEdges=0;
        DMO.TextOnNodes=0;
        DMO.TextOnTriangles=0;
        DMO.TextOnInternalNodes=0;
        DMO.TextOnDirichletNodes=0;
        DMO.TriangleQuality=0;
        DMO.Solution=[];
        DMO.TriContour=0;
        DMO.Offset=0;
        DMO.ShowMesh=1;
        DMO.Periodic=0;
        ForceShowMesh=0;
        n=length(OptionsList);
        for k=1:n
            if isnumeric(OptionsList{k})
                if k~=1
                    error('drawMeshOptions:SolutionMustbeTheFirstParameter',...
                        'The solution vector must be the first parameter');
                end
                DMO.Solution=OptionsList{k};
                DMO.ShowMesh=0;
                ForceShowMesh=0;
            else
                switch lower(OptionsList{k})
                    case {'e', 'edge','edges'}
                        if ~isempty(DMO.Solution)
                            error('drawMeshOptions:IncompatibleOptions',...
                                'You cannot plot the solution AND the numbers of the edges');                            
                        end
                        DMO.TextOnEdges=1;
                    case {'n', 'node', 'nodes'}
                        if ~isempty(DMO.Solution)
                            error('drawMeshOptions:IncompatibleOptions',...
                                'You cannot plot the solution AND the numbers of the nodes');                            
                        end
                        DMO.TextOnNodes=1;
                    case {'t','triangle','triangles'}
                        if ~isempty(DMO.Solution)
                            error('drawMeshOptions:IncompatibleOptions',...
                                'You cannot plot the solution AND the numbers of the triangles');                            
                        end
                        DMO.TextOnTriangles=1;
                     case {'i', 'internalnode','internalnodes'}
                        if ~isempty(DMO.Solution)
                            error('drawMeshOptions:IncompatibleOptions',...
                                'You cannot plot the solution AND the numbers of the internal nodes');                            
                        end
                         DMO.TextOnInternalNodes=1;  

                     case {'d', 'dirichlet','dirichletnodes','dirichletnode'}
                        if ~isempty(DMO.Solution)
                            error('drawMeshOptions:IncompatibleOptions',...
                                'You cannot plot the solution AND the numbers of the Dirichlet nodes');                            
                        end
                         DMO.TextOnDirichletNodes=1;  
                     case {'q', 'trianglequality','quality'}
                        if ~isempty(DMO.Solution)
                            error('drawMeshOptions:IncompatibleOptions',...
                                'You cannot plot the solution AND the quality of the triangles');                            
                        end
                        DMO.TriangleQuality=1;  
                     case {'c', 'tricontour','contour'}
                        DMO.TriContour=1;  
                        DMO.ShowMesh=0;
                        ForceShowMesh=0;
                     case {'o', 'offset'}
                        DMO.Offset=1;  
                    case {'h', 'hidemesh'}
                        DMO.ShowMesh=0;
                        ForceShowMesh=0;
                    case {'s', 'showmesh'}
                        ForceShowMesh=1;
                    case {'p','periodic'}
                        if ~isempty(DMO.Solution)
                            error('drawMeshOptions:IncompatibleOptions',...
                                'You cannot plot the solution AND the numbers of the periodic nodes');                            
                        end                       
                        DMO.Periodic=1;
                    otherwise
                        %kend=k;
                        error('mesh2d:drawMeshOptions:UnknownOption',['Unknown option: ' OptionsList{k}]);
                end
            end
            if ForceShowMesh
                DMO.ShowMesh=1;
            end
        end
    end
    
function [x,y]=GetInternalPoint(X,Y,IsHole)
    xc=(X(1)+X(2))/2;
    yc=(Y(1)+Y(2))/2;
    vx=Y(1)-Y(2);
    vy=X(2)-X(1);       
    if ~IsHole
        vx=-vx;
        vy=-vy;
    end
    %starting value for t: it cannot be too large, to avoid moving outside
    %the region!
    %We cold start with t=eps(min(abs([xc,yc])));   
    %but if xc or yc is zero (or both of course), t is 4.9407e-324, which
    %seems to cause problems to Triangle. So let's use a larger initial
    %value
    t=eps(min(abs([xc,yc])))*256;   
    k=0;
    while k<1024*2,
        x=xc+t*vx;
        y=yc+t*vy;        
        if x~=xc || y~=yc;
            break;
        end
        k=k+1;
        t=t*2;
    end
end

  
    
end
end
