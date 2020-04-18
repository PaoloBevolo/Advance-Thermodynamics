classdef border < matlab.mixin.Copyable
%Define the coordinates of the nodes of the Regions borders and the Boundary conditions associated to each edge

%Version 2017.1
%Copyright 2014-2017 Paolo Bardella
properties (SetAccess=protected)
    X
    Y
    Hole
end
properties
    Bc
end
%%%%%%%%%%%%%%%%%
methods
    function obj=border(varargin)
        %border creates a new border object. This function is generally 
        %called by methods of the shape class.
        if nargin==0
            return 
        elseif nargin==4
            x=varargin{1};x=x(:);
            y=varargin{2};y=y(:);
            if  ~(length(x)==length(y) && isnumeric(x) && isnumeric(y) && isreal(x) && isreal(y))
                error('border:InvalidInputParameter',...
                'The supplied x and y coordinate vectors must have the same length');
            end
            hole=varargin{3};
            if ~(isnumeric (hole) && isscalar(hole) &&( hole==0 || hole==1))
                error('border:InvalidInputParameter',...
                    'The ''hole'' parameter must be 0 or 1');
            end            
            BoundaryConditions=varargin{4};
            if  ~isa(BoundaryConditions,'boundaries')
                error('border:InvalidInputParameter',...
                    'The 4th parameter must be a border object');
            end
            obj.X=x;
            obj.Y=y;
            obj.Hole=hole;
            %ugly, but no other solution, since Boundaries has no public
            %constructor, which would be called by repmat
            obj.Bc=BoundaryConditions(ones(size(x)));
            obj=obj.snap();
        end
    end
    
    function res=isEmpty(obj)
    %isEmpty returns a boolean value indicating if a shape has nodes or not
        res=isempty(obj.X);
    end
    function insertNode(obj,nodeIndex,xy)
    %insertNode insert additional nodes in the border, after a selected one
    %obj=obj.InsertNode(n,xy) inserts after the n-th node of the border and
    %before node n+1, the nodes specified in the two-columns matrix xy,
    %where the first column indicates the x-coordinates and the secondo
    %columns contains the y-coordinates of the new nodes.
    %
    %IMPORTANT: It is useless to call this function simply as
    %bor.insertNode(nodeIndex,xy), since in MATLAB just a copy of bor is modified, 
    %and the original bor is not changed.
    %Instead, always assign the result to bor itsself, or to another
    %variable
    %bor=bor.insertNode(nodeIndex, xy);
    %A warning is generated if the function output is not assigned to a 
    %variable
    %
    %Example
    %   S=shapes.rect;
    %   S.Borders=S.Borders(1).insertNode(1,[0,0]);
        if length(obj)>1
            error('border:InsertNode:InvalidInputParameter',...
                'The first parameter must be a single Region object, not a vector');
        end
        if ~(isnumeric(nodeIndex) && isscalar(nodeIndex) &&isreal(nodeIndex) && nodeIndex>=1 && nodeIndex<=length(obj.X))
            error('border:InsertNode:InvalidInputParameter',...
                'The nodeIndex must be a real scalar number ranging between 1 and the current number of nodes');
        end
        obj.X=[obj.X(1:nodeIndex); xy(:,1); obj.X(nodeIndex+1:end)];
        obj.Y=[obj.Y(1:nodeIndex); xy(:,2); obj.Y(nodeIndex+1:end)];
        obj.Bc=[obj.Bc(1:nodeIndex); obj.Bc(nodeIndex*ones(size(xy,1),1)); obj.Bc(nodeIndex+1:end)];
          
    end    
    function disp(obj)
        obj.info();
    end
    function info(obj, verboseLevel)
        if nargin==1
            verboseLevel=0;
        end
        for p=1:numel(obj)
            str=['    Border ' num2str(p) ': ' nodes2str(length(obj(p).X))];
            if obj(p).Hole==1
                str=[str ', hole'];                                      %#ok
            end
            disp(str);
            if verboseLevel> 2
                disp ([obj(p).X(:), obj(p).Y(:)]);
            end
        end
    end    
    
    function hg = draw(borders,Opt,col,RegionIndex)
    MyNarginchk(nargin,1,4);
    [~,I]=sort([borders(:).Hole]);
        borders=borders(I);
        if nargin<2
            Opt=struct();
        end
        if nargin<3
            col='r';
        end
        if nargin<4
            RegionIndex=1;
        end
         if isfield(Opt,'Colors')
             UseColors=Opt.Colors;
         else
             UseColors=1;
         end

         %XX=zeros(sum(borders.edgesCount())+numel(borders),1);
         %YY=zeros(size(XX));
         %startpos=1;    
         hg=zeros(length(borders),1);
         for p=1:length(borders)
                if numel(borders(p).X)>0
                    if borders(p).Hole
                        hg(p)=patch(borders(p).X, borders(p).Y,'w','HandleVisibility','off','tag','hole');
                    else
                        hg(p)=patch(borders(p).X, borders(p).Y,col,'HandleVisibility','off','tag','solid');
                        if UseColors==0                        
                            set(hg(p),'facecolor','none');
                        end
                    end            
                        
                end
         end          
         

        if isfield(Opt,'Boundary') && Opt.Boundary==1
            Robin=[];
            Dirichlet=[];
            Neumann=[];
            Periodic=[];
            Continuity=[];
            for p=1:length(borders)
                for n=1:numel(borders(p).X)
                    x1=borders(p).X(n);
                    y1=borders(p).Y(n);
                    m=n+1;
                    if m>numel(borders(p).X)
                        m=1;
                    end
                    x2=borders(p).X(m);
                    y2=borders(p).Y(m);
                    switch borders(p).Bc(n).Condition
                        case 'R'
                            Robin=[Robin; x1,y1;x2,y2; NaN NaN];%#ok
                        case 'N'
                            Neumann=[Neumann; x1,y1;x2,y2; NaN NaN]; %#ok
                        case 'D'
                            Dirichlet=[Dirichlet; x1,y1;x2,y2; NaN NaN]; %#ok
                        case 'C'
                            Continuity=[Continuity; x1,y1;x2,y2; NaN NaN]; %#ok
                        case 'P'
                            Periodic=[Periodic; x1,y1;x2,y2; NaN NaN]; %#ok
                    end
                end
            end
            if ~isempty(Dirichlet)
                plot(Dirichlet(:,1),Dirichlet(:,2),'y','linewidth',3,'HandleVisibility','off','tag','dirichlet');
            end
            if ~isempty(Neumann)
                plot(Neumann(:,1),Neumann(:,2),'b--','linewidth',3,'HandleVisibility','off','tag','neumann');
            end
            if ~isempty(Robin)
                plot(Robin(:,1),Robin(:,2),'r-.','linewidth',3,'HandleVisibility','off','tag','robin');
            end
            if ~isempty(Continuity)
                plot(Continuity(:,1),Continuity(:,2),'g-.','linewidth',3,'HandleVisibility','off','tag','continuity');
            end
            if ~isempty(Periodic)
                plot(Periodic(:,1),Periodic(:,2),'m-.','linewidth',3,'HandleVisibility','off','tag','periodic');
            end
        end
        if isfield(Opt,'TextOnNodes') && Opt.TextOnNodes           
            for p=1:length(borders)
                for n=1:numel(borders(p).X)
                    x=borders(p).X(n);
                    y=borders(p).Y(n);
                    if ~isempty(RegionIndex)
                        text(x,y, [num2str(RegionIndex) '\\' num2str(p) '\\' num2str(n)],...
                            'verticalalignment','middle','horizontalalignment','center');
                    elseif length(borders)>1
                        text(x,y, [num2str(p) '\\' num2str(n)],...
                            'verticalalignment','middle','horizontalalignment','center');
                    else
                        text(x,y, num2str(n),...
                            'verticalalignment','middle','horizontalalignment','center');
                    end
                end
            end            
        end

        if isfield(Opt,'TextOnEdges') &&Opt.TextOnEdges     

            for p=1:length(borders)
                for n=1:numel(borders(p).X)
                    x1=borders(p).X(n);
                    y1=borders(p).Y(n);
                    m=n+1;
                    if m>numel(borders(p).X)
                        m=1;
                    end
                    x2=borders(p).X(m);
                    y2=borders(p).Y(m);
                    plot([x1,x2],[y1 y2],'x','HandleVisibility','off','linestyle','none');
                    if ~isempty(RegionIndex)
                        ht=text((x1+x2)/2,(y1+y2)/2, [num2str(RegionIndex) '\\' num2str(p) '\\' num2str(n)]);

                    elseif numel(borders)>1
                        ht=text((x1+x2)/2,(y1+y2)/2, [num2str(p) '\\' num2str(n)]);
                    else
                        ht=text((x1+x2)/2,(y1+y2)/2, num2str(n));
                    end
                    if y1==y2
                        if x1>x2
                            set(ht,'verticalalignment','bottom','horizontalalignment','center','rot',0.1);
                        else
                            set(ht,'verticalalignment','top','horizontalalignment','center','rot',0.1);
                        end
                    elseif x1==x2
                        if y1>y2
                            set(ht,'verticalalignment','bottom','horizontalalignment','center','rot',89);
                        else
                            set(ht,'verticalalignment','bottom','horizontalalignment','center', 'rot',-89);
                        end
                    else
                        hAngle=atan((y2-y1)/( x2-x1))/pi*180;
                        if x1>x2
                            set(ht,'verticalalignment','top','horizontalalignment','center', 'rot',hAngle);
                        else
                            set(ht,'verticalalignment','bottom','horizontalalignment','center', 'rot',hAngle);
                        end
                    end


                end
            end
        end
    end
    function obj=mrdivide(b, val)
        if ~isreal(val)
            error('border:mrdivide:InvalidInputParameter',...
                'Error: the second parameter should be a real scalar, or a 2 componenets real vector');
        end
        if isscalar(val)
           val=[val, val];
        elseif numel(val)~=2
            error('border:mrdivide:InvalidInputParameter',...
                'Error: the second parameter should be a real scalar, or a componenets real vector');
        end
        obj=copy(b);
        for p=1:length(obj)
                obj(p).X=obj(p).X/val(1);
                obj(p).Y=obj(p).Y/val(2);
                obj(p)=snap(obj(p));  
        end 
    end
    function obj=mtimes(b, val)
        if ~isreal(val)
            error('border:mtimes:InvalidInputParameter',...
            'Error: the second parameter should be a real scalar, or a 2 componenets real vector');
        end
        if isscalar(val)
           val=[val, val];
        elseif numel(val)~=2
            error('border:mtimes:InvalidInputParameter',...
                'Error: the second parameter should be a real scalar, or a componenets real vector');
        end
        obj=copy(b);
        for p=1:length(obj)
            obj(p).X=obj(p).X*val(1);
            obj(p).Y=obj(p).Y*val(2);                
            obj(p)=snap(obj(p));  
        end 
        
    end    
    function obj=plus(b1,b2)
        if isnumeric (b2) && isreal(b2)  && numel(b2)==2
            obj=translate(b1,b2);
        elseif isa(b2,'border')
            obj=combine(b1,b2,3);
        else
            error('border:plus:InvalidInputParameter',...
                'Error: to a border object , you can add another Border object or a real two components vector');
        end
    end
    function obj=minus(b1,b2)
        if isnumeric (b2)&& isreal(b2)  && numel(b2)==2
            obj=translate(b1,-b2);
        elseif isa(b2,'border')
            obj=combine(b1,b2,0);
        else
            error('border:minus:InvalidInputParameter',...
                'Error: to a border object, you can subtract another Border object or a real two components vector');
        end
    end
    function obj=and(b1,b2)
        if isa(b2,'border')
            obj=combine(b1,b2,1);
        else
            error('border:and:InvalidInputParameter',...
                'Error: you can calculate the logic intersection between two border objects');
        end
    end   
    
    function obj = rotate(borders,angle)
    %Rotates the points of a shape by the supplied angle
    %Input:
    %   Sh:     Region or vector of shapes to be rotated
    %   angle:  rotation angle in degree
    %Output:
    %   Sh:     The rotated shape(s)
    %REMARKS (very important)
    % The rotation center is always (0,0) and not the shape's centroid.
    MyNarginchk(nargin,2,2);

    if ~(isnumeric(angle) && isreal(angle))
        error('shape:rotate:InvalidInputParameter',...
            'The angle must be a scalar real number');
    end

    %wrap and check angle
    angle=mod(angle, 360);
    if(angle==0)
        return;
    end
    %convert to rad
    angle = angle*pi/180;
    %build rotation matrix
    Rot = [ cos(angle), sin(angle)
         -sin(angle), cos(angle)];
    %loop on the shapes
    obj=borders;
    for p=1:numel(borders)
        node=[borders(p).X.';borders(p).Y.'];
        node = Rot*node;
        obj(p).X=node(1,:).';
        obj(p).Y=node(2,:).';
        obj=snap(obj);   
    end
    
    end
    function A=area(obj)
        A=zeros(size(obj));
        for k=1:numel(obj)
            A(k)=polyarea(obj(k).X,obj(k).Y);
        end
    end    
    function [Mx,My]=min(obj)
        Mx=zeros(size(obj));
        My=zeros(size(obj));
        for p=1:numel(obj)
            Mx(p)=min([obj(p).X]);
            My(p)=min([obj(p).Y]);
        end        
    end    
    function [Mx,My]=max(obj)
        Mx=zeros(size(obj));
        My=zeros(size(obj));
        for p=1:numel(obj)
            Mx(p)=max([obj(p).X]);
            My(p)=max([obj(p).Y]);
        end        
    end    
    
    function Num=edgesCount(obj)
        Num=zeros(size(obj));
        for p=1:numel(obj)
            Num(p)=numel(obj(p).X);
        end        
    end    
    
    function obj=snap(obj,snapFactor)
        if nargin==1
            snapFactor=[grid2D.stepX,grid2D.stepY];
        end
        if ~isempty(snapFactor)
            for kp=1:length(obj)
                if isempty(obj(kp).X)
                    continue;
                end
                x=border.doSnap(obj(kp).X,snapFactor(1));
                y=border.doSnap(obj(kp).Y,snapFactor(2));
                x=[x(:); x(1)];
                y=[y(:); y(1)];
                d2=diff(x).^2+diff(y).^2;
                Keep=d2>=1e-12*max(d2);
                obj(kp).X=x(Keep);
                obj(kp).Y=y(Keep);
                obj(kp).Bc=obj(kp).Bc(Keep);                
                if numel(obj(kp).X)<2
                    obj(kp).X=[];
                    obj(kp).Y=[];
                    obj(kp).Bc=[];
                end
            end
        end
    end
    
    function obj=set.Bc(obj, value)
        if length(value)~=length(obj.X) %#ok there is no other way...
            error('border:bc:WrongInput','The BC vector must have the same length has the coordinates vectors');
        end
        obj.Bc=value;
    end
end

%%%%%%%%%%%%%%%%%
methods (Static)
    function bc = defaultBc()
        bc=boundaries.dirichlet(0);
    end    
    function v=doSnap(v, factor)
        v=round(v/factor)*factor;
    end    
end   
%%%%%%%%%%%%%%%%%
methods (Access=protected)
    function obj=sortBorders(obj)
        [~, index]=sort([obj(:).Hole]);
        obj=obj(index);
    end
    
    function obj=translate(b1,v)   
        obj=copy(b1);
        dx=v(1);
        dy=v(2);
        for k=1:numel(b1)
            obj(k).X=b1(k).X+dx;
            obj(k).Y=b1(k).Y+dy;
        end
        obj=sortBorders(obj);
        obj=snap(obj);        
    end    
    
    function obj=combine(b1,b2,code)                            
        for k=length(b1):-1:1
            b1s(k)=struct('x',b1(k).X,'y',b1(k).Y,'hole',b1(k).Hole);
        end
        for k=length(b2):-1:1
            b2s(k)=struct('x',b2(k).X,'y',b2(k).Y,'hole',b2(k).Hole);
        end        
        tmp=PolygonClip(b1s,b2s,code);
        for k=length(tmp):-1:1
            obj(k)=border(tmp(k).x, tmp(k).y, tmp(k).hole, border.defaultBc);
        end
        obj=sortBorders(obj);
        obj=snap(obj);        
    end
end
end
function str=nodes2str(n)
    if n==1
        str= '1 node';
    else
        str= [num2str(n) ' nodes'];
    end
end