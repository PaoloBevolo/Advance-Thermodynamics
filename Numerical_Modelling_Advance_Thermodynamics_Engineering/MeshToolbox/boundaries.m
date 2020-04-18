classdef boundaries
%Class used to define the boundary conditions (B.C.s) of a shape object

%Version 2017.1
%Copyright 2014-2017 Paolo Bardella
properties (SetAccess=protected)
    Condition
    Value
    Function
end
methods(Static)
    function obj = dirichlet(d)      
    %dirichlet boundary condition        
        if nargin==0
            d=0;
        end        
        val=[];
        func=[];
        if isscalar(d) && isreal(d)
            val=d;
        elseif ~isa(d,'function_handle') || nargin(d)~=2
            error('boundaries:dirichlet:InvalidInputArgument','Dirichlet conditions require a real scalar value or a function of two variables(x,y)')
        else
            func=d;
        end
        obj=boundaries('D',val, func);
    end
    function obj = neumann(d)    
    %neumann boundary condition
        if nargin==0
            d=0;
        end
        val=[];
        func=[];
        if isscalar(d) && isreal(d)
            val=d;
        elseif ~isa(d,'function_handle') || nargin(d)~=2
            error('boundaries:neumann:InvalidInputParameter','Neumann conditions require a real scalar value or a function of two variables(x,y)')
        else
            func=d;
        end
        obj=boundaries('N',val, func);
    end
    
    function obj=periodic(fun)
    %periodic boundary conditions
    %fun is a function handler indicating the periodicity law
        if ~isa(fun,'function_handle') || nargin(fun)~=2
            error('boundaries:periodic:InvalidInputParameter',...
                'Periodic conditions require a function of two variables(x,y) indicating the position of the linked edges')
        end
        obj=boundaries('P',[],fun);
    end      
    function obj=robin(h,g)
    %robin boundary conditions
        if nargin<2
            error('boundaries:robin:InvalidInputParameter','Robin conditions require a couple of real scalar values or functions of two variables(x,y)')
        end        
        val=[];
        func=[];
        if isscalar(h) && isreal(h) && isscalar(g) && isreal(g)
            val=[h,g];
        elseif ~isa(h,'function_handle') || nargin(h)~=2 || ~isa(g,'function_handle') || nargin(g)~=2
            error('boundaries:robin:InvalidInputParameter','Robin conditions require a couple of real scalar values or functions of two variables(x,y)')
        else
            func={h,g};
        end
        obj=boundaries('R',val,func);
    end
    function obj = continuity()  
    %continuity boundary condition
        obj=boundaries('C',[],[]);
    end
end
methods 
    function display(Bc)
        for k=1:length(Bc)
            switch Bc(k).Condition
                case {'D'}
                    str='Dirichlet';
                case {'N'}
                    str='Neumann';
                case {'R'}
                    str='Robin';
                case {'P'}
                    str='Periodic';
                case {'C'}
                    str='Continuity';
                otherwise
                    str='Unknown';
            end
            if ~isempty(Bc(k).Value)
                disp([str ', ' num2str(Bc(k).Value)]);
            elseif ~isempty(Bc(k).Function)
                disp([str ', ' char(Bc(k).Function)]);
            end
        end        
    end
    function v = evaluate( Bc, xy )

if ~isempty(Bc.Value)    
    if isscalar(xy)
        v=Bc.Value{:};
    else
        v=repmat(Bc.Value, size(xy,1),1);
    end
else
    for f=length(e.Fun):-1:1
        v(f,:)=Bc.Fun{f}(xy(:,1),xy(:,2)); 
    end
end
    end
end

methods(Access=protected)
    function obj = boundaries(condition, val, func)
        obj.Condition=condition;
        obj.Value=val;
        obj.Function=func;
    end
end
end