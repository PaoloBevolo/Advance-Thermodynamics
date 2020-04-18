classdef region < matlab.mixin.Copyable
%Class defining a generic region object

%Version 2017.2
%Copyright 2014-2017 Paolo Bardella
properties
    %border object (vector or scalar) describing the region domain
    Borders 
    %symbolic name of the region    
    Name
end
properties (SetAccess=private)    
    %struct of additional material properties assigned using addProperty,
    %removed using removeProperty and evaluated using evaluateProperty
    Properties  
end  
methods
    function Sh=region(varargin)
        %region build a region object, which is a handly way to collect the
        %information regarding 
        % * the geometrical domain of the region
        % * the physical properties of the material of which the region is
        %   made
        %The geometrical properties are described in the Borders class property,
        %which is an object of the border class.
        %
        %The physical properties are described in the Properties class
        %property, which is a structure containing the name of the property
        %and its value (which can be a scalar or a function of x and y)
        %
        %region(x,y) creates a region object described by the points
        %contianed in the x and y vectors.
        %
        %region(sh) creates a copy of the supplied region object sh
        %
        %region(bor) creates region based on the border object bor
        %
        %Additional parameters are possible as name/value pairs.
        %region(...,'name',RegionName) uses the string RegionName as name of the
        %  region, instead of the automatically generated name. In the
        %  string RegionName, the symbol '#' it replaced by the a progressive
        %  number of generated regions.
        %region(...,'boundary',boundaryObject) assignes the boundaries object 
        %  boundaryObject as boundary condition to all the edges of the shap. 
        %  If this parameter is omitted, Homogeneous Dirichet boundary 
        %  conditions are imposed
        %All the other name/value pairs are considered as physical
        %properties and passed to the addProperty function.
        %
        %If area of the created region object is 0, a warning is issued
        %
        %Examples
        % %triangular region
        % Sh=Region([1 2 3],[0 1 0]);        
        % %region object named 'trianglexxx', where xxx is the number of region 
        % %objects already generated
        % Sh=Region([1 2 3],[0 1 0],'name','triangle#'); 
        % %region object with Neumann boundary conditions
        % Sh=Region([1 2 3],[0 1 0],'boundary',boundaries.neumann(0)); 
        % Sh=Region([1 2 3],[0 1 0],'boundary',boundaries.neumann(0)); 
        
        %calculate the number of parameters BEFORE the first string, which
        %is where the name/value pairs start
        if nargin==0
            N=0;
        else
            for N=1:nargin
                if ischar(varargin{N})
                    N=N-1;      %#ok
                    break;
                end
            end
        end
        namefound=false;
        bcsfound=false;
        defautbcs=region.defaultBoundary();
        %search for name/pair values indicating the region name ('name') and
        %default boundary condition ('boundary')
        for k=nargin-1:-2:N+1,
            if ~ischar(varargin{k})
                error('region:region:InvalidPairValueInputs', ...
                    'The optional additional input arguments must be provided as name-value pairs');
            end
            if strcmpi(varargin{k}, 'name') && namefound==false
                %'name' was found for the first time. Since k starts from
                %nargin, this is the last occurrence in the name/value list
                name=varargin{k+1};
                if ~ischar(name)
                     error('region:region:InvalidInputParameter', ...
                        'The name assigned to a region must be a literal string');
                end
                %if it contains '#', replace it by the overall number of
                %already generated regions
                if any(name=='#')
                    name=strrep(name,'#',num2str(region.regionIndex()));
                end
                Sh.Name=name;
                namefound=true;                
            elseif strcmpi(varargin{k}, 'boundary') && bcsfound==false
                %'boundary' was found for the first time
                if ~isa(varargin{k+1},'boundaries')
                     error('region:region:InvalidInputParameter', ...
                        'The default boundary assigned to a region must be a boundaries object');
                end
                bcsfound=true;
                defautbcs=varargin{k+1};
            end
        end
        %look at the non pair/value parameters
        switch N
            %default constructor
            case 0
            case 1
                if isa(varargin{1},'region')
                    Sh=varargin{1};
                elseif isa(varargin{1},'border')
                    Sh.Borders=varargin{1};
                else
                    error('region:region:InvalidInputParameter', ...
                        'It was not possible to extract region informations from the first supplied parameter');
                end
                Sh=Sh.snap();
            case 2
                x=varargin{1};x=x(:);
                y=varargin{2};y=y(:);
                if  length(x)==length(y) && isnumeric(x) && isnumeric(y) && isreal(x) && isreal(y) 
                    Sh.Borders=border(x,y,0,defautbcs);
                else
                    error('region:region:InvalidInputParameter', ...
                        'The supplied x and y coordinate vectors must have the same length');
                end
            otherwise
                error('region:region:InvalidNumberOfInputs', ...
                    'The function expects 1 or 2 input parametrers follwed by optional name-value pairs');
        end
        if N>0
            %if the user, or some region operation, generated any empty border, remove it!
            ToDelete=false(length(Sh.Borders),1);
            for kp=1:length(Sh.Borders)
                if Sh.Borders(kp).isEmpty()
                    ToDelete(kp)=true;
                end
            end
            Sh.Borders(ToDelete)=[];
            %if no border is left, issue a warning
            if isempty(Sh.Borders)
                warning('region:region:Emptyregion', ...
                    'The created region object has no borders');
            end
            %assign the physical properties
            
            for k=N+1:2:length(varargin)
                if k+1>length(varargin)
                     error('region:region:namedparamterrequired', ...
                ['The parameter named ' varargin{k} ' has no corresponding value']);
                end
                if ~strcmpi(varargin{k}, 'name') && ...
                    ~strcmpi(varargin{k}, 'defaultbcs')
                    Sh=Sh.addProperty(varargin{k},varargin{k+1});
                end
            end
            %if the name was not supplied, use the default one
            if isempty(Sh.Name)
                Sh.Name=region.defaultregionName(region.regionIndex());
            end
        end
        
        %if no area, issue a warning
        %%if (isempty(Sh.Borders) && nargout>0) || Sh.area()==0
        if isempty(Sh.Borders)  || Sh.area()==0
            warning('region:region:RegionWithZeroArea', ...
                'The region area is equal to zero');
        end
    end
    
    function disp(Sh)
        Sh.info();
    end
    
    function info(Sh, verboseLevel)
    %display shows textual infos on the provided regions.
    %Sh.display() prints info on the regions and related borders and
    %   properties according to the verbose level specified using the
    %   verbose method
    %Sh.display(verboseLevel) forces a verbose level of printing
    %The greater is the verbose level, the higher is the amount of
    %information printed.
        if nargin==1 | (nargin==2 & ischar(verboseLevel))
            verboseLevel=region.verbose;
        end
        %Count the total number of borders and nodes
        BordersCount=0;
        NodesCount=0;
        for s=1:numel(Sh)
            BordersCount=BordersCount+numel(Sh(s).Borders);
            for p=1:numel(Sh(s).Borders)
                NodesCount=NodesCount+length(Sh(s).Borders(p).X);
            end
        end
        %only one region
        if numel(Sh)==1
            disp([class(Sh) ' object; name: ' Sh.Name '; ' ...
                borders2str(BordersCount) ' and ' nodes2str(NodesCount)]);
        else %vactor of regions
            disp([num2str(size(Sh,1)) 'x' num2str(size(Sh,2)) ' ' class(Sh) ' objects; ' ...
                borders2str(BordersCount) ' and ' nodes2str(NodesCount) ' in total']);
        end

        if verboseLevel>0
            %save current output format, then use compact format to save space
            f=get(0,'FormatSpacing');
            format('compact'); 
            for s=1:numel(Sh)
                disp(['region ' num2str(s) ' (' Sh(s).Name '): ' borders2str(numel(Sh(s).Borders)) ]);
                if ~isempty(Sh(s).Properties) && ~isempty(fieldnames(Sh(s).Properties))
                    disp(Sh(s).Properties);
                end
                if verboseLevel>1
                    info(Sh(s).Borders, verboseLevel);
                end
            end
            set(0,'FormatSpacing',f); %restore previous output format
        end
    end
    
    function draw(Sh, varargin)
    %draw draws the region objects using a 2D graph
    %Different behaviours can be obtained according to the optional
    %string parameters which are specified as draw('param1','param2',...)
    %
    %Each string can be 
    %   e, edge, edges:     to print the number of each edge
    %   n, node, nodes:     to print the number of each node
    %   bw:                 to display only the regions borders
    %                       without filling
    %   bc:                 to display the boundary conditions. The
    %                       following color/line style code is used:
    %                       Dirichlet BCs: yellow line
    %                       Neumann BCs: blue dashed line
    %                       Robin BCs: red dashed-dotted line
    %                       Periodic BCs: green dashed-dotted line 
    %                       Continuity BCs: magenta dashed-dotted line
    % and they can be combined in any sequence.
    %
    %Examples
    %Sh=regions.rect();
    %Sh.draw('e');
        Opt=region.drawOptions(varargin);
        %save the hold condition of the current axis to restore it at the
        %end
        currentstate=get(gca,'nextplot');
        hold('on');
        %hobj=zeros(numel(Sh),1);
        col=lines(numel(Sh));
        %for each region object
        for s=1:numel(Sh)            
            borders=Sh(s).Borders;
            if numel(Sh)>1
                RegionIndex=s;
            else
                RegionIndex=[];
            end
            hobj=borders.draw(Opt, col(s,:), RegionIndex);           
            htmp=findall(hobj,'tag','solid');
            set(htmp(1),'DisplayName',Sh(s).Name,'handleVisibility', 'on');
        end 
        h=findall(gca,'tag','dirichlet');
        if ~isempty(h)
            set(h(1),'HandleVisibility','on','DisplayName','Dirichlet');
        end
        h=findall(gca,'tag','neumann');
        if ~isempty(h)
            set(h(1),'HandleVisibility','on','DisplayName','Neumann');
        end
        h=findall(gca,'tag','robin');
        if ~isempty(h)
            set(h(1),'HandleVisibility','on','DisplayName','Robin');
        end
        h=findall(gca,'tag','periodic');
        if ~isempty(h)
            set(h(1),'HandleVisibility','on','DisplayName','Periodic');
        end
        h=findall(gca,'tag','continuity');
        if ~isempty(h)
            set(h(1),'HandleVisibility','on','DisplayName','Continuity');
        end   
      
        %sort the graphic objects: in the background the patches, then the line, finally the
        %text
        patches=findall(gca,'type','patch');
        for k=1:length(patches),
            try
            uistack(patches(k),'bottom');
            catch
            end
        end

            
        patches=findall(gca,'type','text');
        
        for k=1:length(patches),
            try            
                uistack(patches(k),'top');
            catch
            end

        end
      
        legend('show');
        set(gca,'nextplot',currentstate);
        [Mx, My]=max(Sh);
        [mx, my]=min(Sh);
        X=xlim();
        Y=ylim();
        a=[X(1) X(2) Y(1) Y(2)];
        StepX=(X(2)-X(1))/20;
        StepY=(Y(2)-Y(1))/20;
        if min(mx)==X(1)
            a(1)=a(1)-StepX;
        end
        if min(my)==Y(1)
            a(3)=a(3)-StepY;
        end
        if max(Mx)==X(2)
            a(2)=a(2)+StepX;
        end
        if max(My)==Y(2)
            a(4)=a(4)+StepY;
        end
        axis(a);
        if isempty(get(get(gca, 'xlabel'),'string'))
            xlabel('x dir [m]');
        end
        if isempty(get(get(gca, 'ylabel'),'string'))
            ylabel('y dir [m]');
        end
        set(gca,'box','on');
    end
    
    function Sh=addProperty(Sh,name, value)
    %addProperty add a physical property to the specified regions
    %Sh=Sh.addProperty('propertyName',propertyValue) adds the physical 
    %   property propertyName to Sh. The assigned value is propertyValue,
    %   which can be either a real vector/matrix, or a function receiving
    %   in input x and y coordinates and returning a real vector/matrix.
    %propertyName must be a valid name for a field in a MATLAB structure
    %
    %If a property called propertyName already exists, it is overwritten
    %without warnings
    %
    %To evaluate the property, use the member function evaluateProperty
    %
    %To remove a property, use the member function removeProperty
    %
    %Example
    % Sh.addProperty('density',1);    
        MyNarginchk(nargin,3,3);
        if ~isvarname(name)
            error('region:addProperty:InvalidProperty', ...
                'The property name must be a valid MATLAB identifier');
        end
        if (isnumeric(value) && isreal(value))
            for k=1:numel(Sh)
                Sh(k).Properties.(name)=value;    
            end
        elseif (isa(value,'function_handle') )
            if nargin(value)~=2
                error('region:addProperty:InvalidProperty', ...
                      'The function handle must have two input arguments');
            end
            for k=1:numel(Sh)
                Sh(k).Properties.(name)=value;
            end
        else
            error('region:addProperty:InvalidProperty', ...
                'The property value must either a real number or a function handle');
        end       
    end
    
    function Sh=removeProperty(Sh,name)
    %removeProperty removes from a region a previously added physical property 
    %Sh=Sh.removeProperty('propertyName') removes the physical property 
    %   propertyName from Sh. 
    %
    %If a property called propertyName does not exists, a warning is
    %issued.
    %
    %IMPORTANT: It is useless to call this function simply as
    %Sh.addProperty('propertyName',propertyValue), since in MATLAB just a 
    %copy of Sh is modified, and the original Sh is not changed.
    %Instead, always assign the result to Sh itsself, or to another
    %variable
    %Sh=Sh.addProperty('propertyName',propertyValue);
    %A warning is generated if the function output is not assigned to a 
    %variable
    %
    %Example
    % Sh=Sh.removeProperty('density');
        MyNarginchk(nargin,2,2);        
        if ~ischar(name) %I simply check for name to be a string.
            error('region:removeProperty:InvalidProperty', ...
                'The property name must be a string');
        end        
        found=false;
        for k=1:numel(Sh)
            if isfield(Sh(k).Properties, name)
                Sh(k).Properties=rmfield(Sh(k).Properties,name);
                found=true; %found at least once in the regions vector
            end
        end
        if ~found %never found
            warning('region:removeProperty:InvalidProperty', ...
                ['The property ' name ' does not exist in the provided region(s)']);
        elseif nargout==0
            warning('region:removeProperty:ResultNotSaved', ...
                'The property has been set, but the result was not saved in a variable');
        end
    end
    
    function res = evaluateProperty(Sh,propertyName, x, y)
    %evaluateProperty evaluates the specified material property 
    %Sh.evaluateProperty(propertyName, xy) evaluate the specified
    %property in the point indicated by the components x and y
    %
    %If a property called propertyName does not exists, an error is
    %issued.        
    %Example
    % Sh=Sh.evaluateProperty('density',0,0);
    
        %additional tests on the property.
        try %faster than isfield
            res=Sh.Properties.(propertyName);    
        catch
             error('region:evaluateProperty:InvalidProperty', ...
                ['The property ' propertyName ' does not exist in the provided region(s)']);
        end
        try
            if ~isnumeric(res)
                res=res(x,y);
            else
                res=repmat(res,size(x,1),1);
            end
        catch Exception                
            error('region:evaluateProperty:InvalidProperty', ...
                [ 'The evaluation of the property ' propertyName ' returned an error:\n' Exception.message]);
        end        
    end
    
    function A=area(Sh)
    %area calculates the ares of a region or a vector of regions.
    %Example
    %Sh=regions.rect();
    %disp(Sh.area)
        A=zeros(size(Sh));
        %for each region object
        for s=1:numel(Sh)
            Area=Sh(s).Borders.area;
            for p=1:length(Sh(s).Borders)
                %add if it's an outern regions
                if Sh(s).Borders(p).Hole==0
                    A(s)=A(s)+Area(p);
                else
                    %subtract if it's an hole
                    A(s)=A(s)-Area(p);
                end            
            end
        end
    end    
    
    function obj=snap(Sh)
    %Snap adapts the regions' borders to the grid2D grid and returns the 
    %corrected regions. 
        obj=Sh;
        for s=1:numel(Sh)
            obj(s).Borders=Sh(s).Borders.snap();
        end
        if nargout==0
            warning('region:snap:ResultNotSaved', ...
                'The region has been snapped, but the result was not saved in a variable');
        end        
    end
       
    function [mx,my]=min(Sh)
    %min returns the minimum of the x and y coordinates of the rectangle
    %with sides parallel to the axes and enclosing the provided region.
    %[mx,my]=min(Sh) when Sh is a vector of regions, returns two column vectors
    %containinig the minimum x and y of each region.
        mx=zeros(size(Sh));
        my=zeros(size(Sh));
        for s=1:length(Sh)
            [mx,my]=Sh(s).Borders.min();
        end
    end
    
    function [Mx,My]=max(Sh)
    %max returns the maximum of the x and y coordinates of the rectangle
    %with sides parallel to the axes and enclosing the provided region.
    %[mx,my]=max(Sh) when Sh is a vector of regions, returns two column vectors
    %containinig the maximum x and y of each region.        
        Mx=zeros(size(Sh));
        My=zeros(size(Sh));
        for s=1:length(Sh)
            [Mx,My]=Sh(s).Borders.max();
        end
    end
    
    function obj=plus(obj1,obj2)
    %plus either adds two regions together, or offsets a region by a
    %specified amount
    %Shout=plus(Sh1, Sh2), where Sh1 and Sh2 are both region objects, 
    %   returns a region (or a vector of regions) which is the union ("sum") 
    %   of the 2D domains described by Sh1 and Sh2. 
    %Shout=plus(Sh1, [x,y]) translate the domain described by the region Sh1
    %   by x in the x direction and y in the y direction.
    %Example:
    %   Sh1=regions.rect();
    %   Sh2=Sh1+[0.5,0.5]
    %   Sh3=Sh1+Sh2;
    %   Sh3.draw();
    
        if isa(obj2,'region')
            if any(size(obj1)~=size(obj2))
                error('region:InvalidInputSize',...
                    ['Input arguments sizes are not consistent: region1 is ' num2str(size(obj1)) '  while region is' num2str(size(obj2))]);
            end
            obj=obj1;
            for s=1:numel(obj1)
                obj(s)=region(obj1(s).Borders+obj2(s).Borders,'name',obj1(s).Name);
                obj(s).Properties=obj1(s).Properties;
            end
        elseif isnumeric (obj2)            
            if numel(obj2)~=2
                error('region:InvalidInputSize',...
                    'The offset must be a 2-elements vector');
            end            
            obj=obj1;
            for s=1:numel(obj1)
                obj(s)=region(obj1(s).Borders+obj2,'name',obj1(s).Name);
                obj(s).Properties=obj1(s).Properties;
            end
        else
            error('region:InvalidInputSize',...
                'The second term must be either a region object or a 2-elements vector');
        end
    end
    
    function obj=minus(obj1,obj2)
    %minus either "subtracts" two regions, or offsets a region by a
    %   specified amount
    %Shout=minus(Sh1, Sh2), where Sh1 and Sh2 are both region objects, 
    %   returns a region (or a vector of regions) which describes the
    %   "relative complement" of Sh2 in Sh1 
    %Shout=minus(Sh1, [x,y]) translate the domain described by the region Sh1
    %   by -x in the x direction and -y in the y direction.
    %Example:
    %   Sh1=regions.rect();
    %   Sh2=Sh1+[0.5,0.5]
    %   Sh3=Sh1-Sh2;
    %   Sh3.draw();        
        if isa(obj2,'region')
            if any(size(obj1)~=size(obj2))
                error('region:InvalidInputSize',...
                    ['Input arguments sizes are not consistent: region1 is ' num2str(size(obj1)) '  while region is' num2str(size(obj2))]);
            end
            obj=obj1;
            for s=1:numel(obj1)
                obj(s)=region(obj1(s).Borders-obj2(s).Borders,'name',obj1(s).Name);
                obj(s).Properties=obj1(s).Properties;
            end
        elseif isnumeric (obj2)            
            if numel(obj2)~=2
                error('region:InvalidInputSize',...
                'The offset must be a 2-elements vector');
            end            
            obj=obj1;
            for s=1:numel(obj1)
                obj(s)=region(obj1(s).Borders-obj2,'name',obj1(s).Name);
                obj(s)=obj(s).snap();
                obj(s).Properties=obj1(s).Properties;
            end
        else
            error('region:InvalidInputSize',...
                'The second term must be either a region object or a 2-elements vector');
        end
    end
    
    function obj=and(obj1,obj2)
    %and returns the intersection of two regions
    %Shout=and(Sh1, Sh2), where Sh1 and Sh2 are both region objects, 
    %   returns a region (or a vector of regions) which describes the
    %   intersection of the 2D domains described by Sh1 and Sh2
    %Example:
    %   Sh1=regions.rect();
    %   Sh2=Sh1+[0.5,0.5]
    %   Sh3=Sh1&Sh2;
    %   Sh3.draw();            
        if isa(obj2,'region')
            if any(size(obj1)~=size(obj2))
                error('region:InvalidInputSize',...
                    ['Input arguments sizes are not consistent: region1 is ' num2str(size(obj1)) '  while region is' num2str(size(obj2))]);
            end
            obj=obj1;
            for s=1:numel(obj1)
                obj(s)=region(obj1(s).Borders&obj2(s).Borders,'name',obj1(s).Name);
                obj(s)=obj(s).snap();
                obj(s).Properties=obj1(s).Properties;
            end
        
        else
            error('region:InvalidInputSize',...
                'The second term must be either a region object or a 2-elements vector');
        end
    end
    
    function obj=mrdivide(Sh, val)
    %mrdivide divides the x and y coordinates of each vertex of the
    %domains described by Sh, by the amount provided in val
    %Shout=Sh/[valx, valy], where valx and valy are both scalar quantities, 
    %   divides the x-coordinates by valx and the y-coordinates by valy
    %Shout=Sh/val, where val is a scalar quantity, 
    %   is equivalen to Shout=Sh/[val, val]
    %Example:
    %   Sh1=regions.rect();
    %   Sh2=Sh1/20;
    %   Sh2.draw();      
        if ~isreal(val)
            error('region:mrdivide:InvalidOperands',...
                'The second parameter should be a real scalar, or a 2 componenets real vector');
        end
        if isscalar(val)
           val=[val, val];
        elseif numel(val)~=2
            error('region:mrdivide:InvalidOperands',...
                'The second parameter should be a real scalar, or a 2 componenets real vector');
        end
        obj=copy(Sh);
        for s=1:numel(Sh)
            obj(s).Borders=obj(s).Borders/val;
            obj(s).Properties=Sh(s).Properties;
        end 
    end
    
    function obj=mtimes(Sh, val)
    %mtimes multiplies the x and y coordinates of each vertex of the
    %domains described by Sh, by the amount provided in val
    %Shout=Sh*[valx, valy], where valx and valy are both scalar quantities, 
    %   multiplies the x-coordinates by valx and the y-coordinates by valy
    %Shout=Shal, where val is a scalar quantity, 
    %   is equivalen to Shout=Sh/[val, val]
    %Example:
    %   Sh1=regions.rect();
    %   Sh2=Sh1*20;
    %   Sh2.draw();      
        
        if ~isreal(val)
            error('region:mtimes:InvalidOperands',...
                'The second parameter should be a real scalar, or a 2 componenets real vector');
        end
        if isscalar(val)
           val=[val, val];
        elseif numel(val)~=2
            error('region:mtimes:InvalidOperands',...
                'The second parameter should be a real scalar, or a 2 componenets real vector');
        end
        obj=copy(Sh);
        for s=1:numel(Sh)
            obj(s).Borders=obj(s).Borders*val;
        end 
    end
    
    function obj = rotate(Sh,angle)
    %Rotates the points of a region by the supplied angle
    %Input:
    %   Sh:     Region or vector of regions to be rotated
    %   angle:  rotation angle in degree
    %Output:
    %   Sh:     The rotated region(s)
    %REMARKS (very important)
    % The rotation center is always (0,0) and not the region's centroid.
    MyNarginchk(nargin,2,2);

    if ~(isscalar(angle) &&isnumeric(angle) && isreal(angle))
        error('region:rotate:InvalidInputParameter',...
            'The angle must be a scalar real number');
    end
    %loop on the regions
    obj=copy(Sh);
    for s=1:numel(Sh)
        obj(s).Borders=obj(s).Borders.rotate(angle);
        obj(s)=obj(s).snap();
    end
    
    end
end
methods (Static)
    function Old=verbose(newLevel)
    %verbose allows to tune the quantity of info printed by the display
    %function.
    %verbose(newLevel) sets the verose level to the scalar non-negative
    %level newLevel. 0 corresponds to the level with less information printed 
    %verbose() returns the current level
    %
    %Example
    %   verbose(1);
    %   current=verbose();
        persistent VerboseLevel;
        if nargin==0
            Old=VerboseLevel;
        elseif isnumeric(newLevel) && isreal(newLevel) && newLevel>=0 && fix(newLevel)==newLevel
            Old=VerboseLevel;
            VerboseLevel=newLevel;
        else
            error('region:InvalidVerboseLevel',' The verbose level must be a non negative integer');
        end
        if isempty(Old)
            Old=0;
        end
    end
    
    function Index=regionIndex(New)
        persistent regionIndex;
        if nargin
            if isnumeric(New) && New>=0 && fix(New)==New
                regionIndex=New;
                Index=New;
            else
                error('region:InvalidregionIndex',' The region index must be a non negative integer');
            end
        else
            if isempty(regionIndex)
                regionIndex=1;
            end
            Index=regionIndex;
            regionIndex=regionIndex+1;
        end
    end
end
methods (Static,Access=protected)
    function DMO=drawOptions(parlist)
        %Create a struct with the drawing parameters for Draw
        %   Available parameters are:
        %   e, edge, edges:       print the number of each edge
        %   n, node, nodes:       print the number of each node
        %   c, color, colors:     display each region with a different random color
        %   bw:         display only region borders
        %   bc:         
        DMO=struct();
        n=length(parlist);
        for k=1:n
            switch lower(parlist{k})
                case {'e', 'edge','edges'}
                    DMO.TextOnEdges=1;
                case {'n', 'node', 'nodes'}
                    DMO.TextOnNodes=1;
                case {'c','color','colors'}
                    DMO.Colors=1;
                case {'bw'}
                    DMO.Colors=0;
                case {'bc'}
                    DMO.Boundary=1;
                otherwise
                    error('region:draw:UnknownOption',['Undefined option: ' parlist{k}]);
            end
        end
    end
    
    function Name=defaultregionName(k)
    %defaultregionName returns the default name of a newly created region.
    %The name is used when no other names are assigned by the user
    %region.defaultName(n) returns 'Sn' if n is a number or a string,
    %otherwise simply 'S'.
        if isnumeric(k)
            Name=sprintf('S%d',k);
        elseif ischar(k)
            Name=['S',k];
        else
            Name='S';
        end
    end
    
    function b=defaultBoundary()    
    %defaultBoundary returns the default Boundary COndition used when a 
    %new region object is created, i.e., an Homogeneous Dirichlet Boundary object 
        b=boundaries.dirichlet(0);
    end
end
end

function str=borders2str(n)
%borders2str converts the number of borders to a string
    if n==1
        str= '1 border';
    else
        str= [num2str(n) ' borders'];
    end
end

function str=nodes2str(n)
%nodes2str converts the number of nodes to a string
    if n==1
        str= '1 node';
    else
        str= [num2str(n) ' nodes'];
    end
end