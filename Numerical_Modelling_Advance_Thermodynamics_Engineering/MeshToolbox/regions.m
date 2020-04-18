classdef regions
%regions simplifies the definition of common geometric regions 

%Version 2017.1
%Copyright 2014-2017 Paolo Bardella
methods (Static, Access=public)
    function R=rect(varargin)
    %rect creates a rectangular region object given its center and dimensions
    %
    %rect() generates a square centered in the origin with unitary side
    %
    %rect(xy) generates a square with unitary side, centered in the position 
    %specified in the 2 components vector xy
    %
    %rect(xy,wh) generates a rectangle, centered in the position specified 
    %in the 2 components vector xy, whose height and width are specified in
    %the 2 components vector wh
    %
    %It is possible to provide additional properties as name/value pairs.
    %
    %rect returns the created region object
    %
    %Examples
    %   Sh=regions.rect();
    %   Sh=regions.rect([1,3]) %square centered in x=1, y=3 with unitary side
    %   Sh=regions.rect([1,3],[2,4] %rectangle centered in x=1, y=3 with heigth=2 and width=4
    
        [N, center, dim] =regions.GetFixedParams({[0,0],1},varargin{:});
        if ~(numel(dim)<3 && isreal(dim) && isnumeric(dim))
            error('regions:rect:InvalidInputParameter','The Dimension must be a real scalar or a vector with 2 components (w and h)');
        end
        if isscalar(dim)
            dim=[dim dim];
        end
        if ~(numel(center)==2 && isreal(center) && isnumeric(center))
            error('regions:rect:InvalidInputParameter','The Center must be a vector with 2 components (x and y)');
        end
        xc=center(1);
        yc=center(2);
        xdim=dim(1)/2;
        ydim=dim(2)/2;
        x=[xc-xdim xc-xdim xc+xdim xc+xdim];
        y=[yc-ydim yc+ydim yc+ydim yc-ydim];
        R = region(x,y,'name','R#', varargin{N+1:end});
    end
    
    function R = rectN(varargin)
    %rectN creates a rectangular region object given two opposites points
    %
    %rectN(xy1,xy2) generates a rectangle have a vertices in the points 
    %indicated by the 2 components vectors xy1 and xy2.
    %
    %It is possible to provide additional properties as name/value pairs.
    %
    %rectN returns the created region object
    %
    %Example
    %   Sh=regions.rectN([1,2],[3,4]);
    
        [N, xy1, xy2] =regions.GetFixedParams({[-0.5,-0.5],[0.5,0.5]},varargin{:});
        if N~=2
            error('regions:rectN:InvalidInputParametersNumber',...
                'This function requires 2 numeric parameters');
        end
        if ~(numel(xy1)==2 && isnumeric(xy1) && isreal(xy1))
            error('regions:rectN:InvalidInputParameter',...
                'The first parameter must be a compenent vector of real numbers (x and y)');
        end
        if ~(numel(xy2)==2 && isnumeric(xy2) && isreal(xy2))
            error('regions:rectN:InvalidInputParameter',...
                'The second parameter must be a compenent vector of real numbers (x and y)');
        end
        x=[xy1(:,1) xy1(:,1) xy2(:,1) xy2(:,1)];
        y=[xy1(:,2) xy2(:,2) xy2(:,2) xy1(:,2)];
        R=region(x,y,'name','R#', varargin{N+1:end});
    end
    
    function C = circle(varargin)
    %circle creates a region whose border is a regular polygon inscripted in 
    %a circle or ellipsis.
    %
    %circle() generates a "circular" region centered in the origin, with unitary
    %diameter. 32 points are used to discretize the circle.
    %
    %circle(xy) generates a "circular" region with unitary diameter, centered 
    %in the position specified in the 2 components vector xy. 32 points are 
    %used to discretize the circle.
    %
    %circle(xy,r) generates a "circular" region, centered in the position specified 
    %in the 2 components vector xy, whose radius is specified by the second 
    %parameter
    %
    %circle(xy,axesxy) generates an "elliptical" region, centered in the position 
    %specified in the 2 components vector xy, whose x and y axes are
    %specified in 2 components vector axesxy. 32 points are used to
    %discretize the ellipsis
    %
    %circle(xy,axesxy,nPoints) uses the specified number of points to
    %discretize the ellipsis.    
    %
    %It is possible to provide additional properties as name/value pairs.
    %
    %circle returns the created region object
    %
    %Examples
    %   Sh=regions.circle();
    %   Sh=regions.circle([1,3],[2,5]) 
    %   Sh=regions.circle([1,3],[2,5],128) 
    
        [N,center,radius,nPoints]=regions.GetFixedParams({[0,0],0.5,32},varargin{:});
        if ~(numel(center)==2 && isnumeric(center) && isreal(center))
            error('regions:circle:InvalidInputParameter',...
                'The first parameter must be a vector with 2 components (x and y)');
        end
        if isscalar(radius)
            radius=[radius,radius];
        end
        if ~(numel(radius)==2 && isnumeric(radius) && isreal(radius) && all(radius>0))
            error('regions:circle:InvalidInputParameter',...
                'The second  parameter must be a positive scalar or a 2-components vector of real numbers');
        end
        if ~(numel(nPoints)==1 && isnumeric(nPoints) && isreal(nPoints) && nPoints>2)
            error('regions:circle:InvalidInputParameter',...
                'The third parameter must be real scalar >2');
        end

        teta=linspace(0,2*pi,nPoints+1);
        %remove the last point since we don't want a closed polygon,
        %otherwise mesh wouldn't work
        teta=teta(1:end-1);
        x=cos(teta)*radius(1)+center(1);
        y=sin(teta)*radius(2)+center(2);
        C=region(x,y,'name','C#', varargin{N+1:end});
    end   
    
    function P = triangle(varargin)
    %Generates a equilateral triangular region
    %Input:
    %   center: 2 components vector indicating the x and y coordinates of the center
    %   side:   scalar indicating the side of the equilater triangle
    %Output:
    %   P:      Object of region class
    [~ , center, side]=regions.GetFixedParams({[0,0],1},varargin{:});
    if numel(center)~=2
        error('Center must be a vector with 2 components (x and y)');
    end
    Sh=regions.circle([0 0],side/sqrt(3),3);
    P=Sh.rotate(-90)+center;
    %P=Rotate(Circle([0 0],side/sqrt(3),3),-90)+center;
    end
    
    function S = sector( varargin )
    %sector creates a region which a sector of elliptical annulus.
    %
    %sector() generates an annulus centered in the origin, with external
    %radius=0.5 and internal radius=0.25. 32 points are used to discretize 
    %the borders.
    %
    %sector(xy) generates a annular region with external radius=0.5 and 
    %internal radius=0.25, centered in the position specified in the 2 
    %components vector xy. 32 points are used to discretize the borders.
    %
    %sector(xy,axesxy) allows to specify the external and internal radii in
    %the two-components vector axesyy
    %
    %sector(xy,axesxy,angles) allows to specify limit the region to the
    %region between the two angles specified (in degrees) in the "angles" vector 
    %
    %sector(xy,axesxy,angles, nPoints) allows to control the number of points 
    %used to approximate the curved domain
    %
    %It is possible to provide additional properties as name/value pairs.
    %
    %sector(xy,axesxy,nPoints) uses the specified number of points to
    %discretize the ellipsis.    
    %
    %sector returns the created region object
    %
    %Examples
    %   Sh=regions.sector();
    %   Sh=regions.sector([1,3],5,2) 
    %   Sh=regions.sector([1,3],5,2,[0, 180],128)     
    %   Sh=regions.sector([1,3],5,2,[0, 90],128,'name','region1','density',1)     
    
        [N,center,radiusExt,radiusInt,angles,nPoints]=...
            regions.GetFixedParams({[0,0],0.5,0.35,[0,360], 32},varargin{:});
        if ~(numel(center)==2 && isnumeric(center) && isreal(center))
            error('regions:sector:InvalidInputParameter',...
                'The first parameter (center) must be a vector with 2 components (x and y)');
        end
        if isscalar(radiusExt)
            radiusExt=[radiusExt,radiusExt];
        end
        if ~(numel(radiusExt)==2 && isnumeric(radiusExt) && isreal(radiusExt) && all(radiusExt)>0)
            error('regions:sector:InvalidInputParameter',...
                'The second  parameter (external radius) must be a positive scalar or a 2-components vector of real numbers');
        end
        if isscalar(radiusInt)
            radiusInt=[radiusInt,radiusInt];
        end
        if ~(numel(radiusInt)==2 && isnumeric(radiusInt) && isreal(radiusInt) && all(radiusInt)>=0)
            error('regions:sector:InvalidInputParameter',...
                'The third  parameter (internal radius) must be a positive scalar or a 2-components vector of real numbers');
        end
        
        if ~(numel(angles)==2 && isnumeric(angles) && isreal(angles))
            error('regions:sector:InvalidInputParameter',...
                'The fourth  parameter (angles) must be a scalar or a 2-components vector of real numbers');
        end
        
        if ~(numel(nPoints)==1 && isnumeric(nPoints) && isreal(nPoints) && nPoints>2)
            error('regions:sector:InvalidInputParameter',...
                'The fifth parameter (number of points) must be real scalar >2');
        end        
        
        if angles(1)==angles(2) ||(angles(1)==0 && angles(2)==360)
            S=regions.circle(center,radiusExt,nPoints,'name','dummy',varargin{N+1:end})-...
                regions.circle(center,radiusInt,nPoints,'name','dummy');
            S.Name=['E' num2str(region.regionIndex)];
        else    
            teta=linspace(angles(1),angles(2),nPoints+1)/360*2*pi;
            xext=cos(teta)*radiusExt(1)+center(1);
            yext=sin(teta)*radiusExt(2)+center(2);
            if radiusInt>0
                teta=teta(end:-1:1);
                xint=cos(teta)*radiusInt(1)+center(1);
                yint=sin(teta)*radiusInt(2)+center(2);
            else
                xint=center(1);
                yint=center(2);
            end
            S=region([xext xint], [yext yint],'name','E#', varargin{N+1:end});
        end
    end
end

methods(Static, Access=protected)
    function [N,varargout]=GetFixedParams(DefaultValues, varargin)
    %Helper function to process user input and assign default values to
    %missing arguments
        if length(DefaultValues)~=nargout-1
            error('regions:GetFixedParams:InvalidInputParameter',...
                'The number of defaults values must equal the number of output parameters');
        end
        N=length(varargin);        
        varargout=DefaultValues;
        for k=1:N
            if ~ischar(varargin{k})
                %'regular' parameter: I copy the value in varargout
                if nargout-1<k
                    error('regions:GetFixedParams:TooManyInputParameters',...
                        'Too many input parameters');
                end
                if isempty(varargin{k})
                    varargout{k}=DefaultValues{k};
                else
                    varargout{k}=varargin{k};
                end                
            else %string: no more 'regular' parameters: use dafault values
                for kk=k:nargout-1
                    varargout{kk}=DefaultValues{kk};
                end
                N=k-1;
                break            
            end
        end    
    end
end
end