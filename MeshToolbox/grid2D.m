classdef grid2D
%Class used to define the the grid2D used to anchor objects

%Version 2017.1
%Copyright 2014-2017 Paolo Bardella
methods (Static, Access=protected)       
    function Step=DefaultStep()
        Step=1e-10;
    end
    function old=IsInitialized(newStatus)
        persistent InitializationStatus;
        if nargin==0
            old=InitializationStatus;
        else
            old=InitializationStatus;
            InitializationStatus=newStatus;
        end            
    end
end
methods (Static)       
    function setSteps(newX, newY)
    %setSteps sets the x and y directions grid2D steps
    %
    %setSteps(step) uses the same step for the x and y directions. If step
    %is an empty vector, the grid2D is removed.    
    %
    %setSteps(stepX, stepY) allows to specify the x and y step separately
    %    
    %Any change affects only the Region objects which are created after
    %calling this function. Existing Region objects are not modified.
    %
    %Example
    %grid2D2D.setStep(1e-4,1e-5);   
    %grid2D2D.setStep([]);       
        if nargin==1
            newY=newX;
        end
        grid2D.stepX(newX);
        grid2D.stepY(newY);
    end
    function old=stepX(new)
    %stepX gets/sets the x direction grid2D step
    %
    %stepX(step) uses the provided input value as step for the x direction. 
    %If step is an empty vector, the grid2D is removed in the x direction.    
    %
    %stepX returns the current grid2D step in the x direction
    %
    %Any change affects only the Region objects which are created after
    %calling this function. Existing Region objects are not modified.
    %
    %Example
    %step=grid2D.stepX;
    %grid2D.stepX(1e-4,1e-5);   
    %grid2D.stepX([]);           
    
        persistent StepXValue;
        if grid2D.IsInitialized
            old=StepXValue;
        else
            old=grid2D.DefaultStep();
        end
        if nargin>0
            if isempty(new) || (isnumeric(new) && isscalar(new) && isreal(new) && new>0)
                old=StepXValue;
                StepXValue=new;
                grid2D.IsInitialized(true);
            else
                error('grid2D:InvalidSnapFactor',' The snap factor must be a non negative number');
            end            
        end
    end
    function old=stepY(new)
    %stepY gets/sets the y direction grid2D step
    %
    %stepY(step) uses the provided input value as step for the y direction. 
    %If step is an empty vector, the grid2D is removed in the y direction.    
    %
    %stepY returns the current grid2D step in the y direction
    %
    %Any change affects only the Region objects which are created after
    %calling this function. Existing Region objects are not modified.
    %
    %Example
    %step=grid2D.stepY;
    %grid2D.stepY(1e-4,1e-5);   
    %grid2D.stepY([]);            
        persistent StepYValue;
        if grid2D.IsInitialized()
            old=StepYValue;
        else
            old=grid2D.DefaultStep();
        end
        if nargin>0
            if isempty(new) || (isnumeric(new) && isscalar(new) && isreal(new) &&new>0)        
                StepYValue=new;
                grid2D.IsInitialized(true);

            else
                error('grid2D:InvalidSnapFactor',' The snap factor must be a non negative number');
            end
        end
    end
    function reset()
    %reset the grid to the original steps
    %
    %Any change affects only the Region objects which are created after
    %calling this function. Existing Region objects are not modified.
    %
    %Example
    %grid2D.reset;
        grid2D.stepX([]);
        grid2D.stepY([]);
        grid2D.IsInitialized(false);
    end 

    function remove()
    %remove the grid in the x and y directions
    %
    %Any change affects only the Region objects which are created after
    %calling this function. Existing Region objects are not modified.
    %
    %Example
    %grid2D.remove;
        grid2D.stepX([]);
        grid2D.stepY([]);
    end 
end
end
    