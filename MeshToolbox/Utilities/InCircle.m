function res=InCircle(xy,xyc,r)
%Returns if a point (x,y) lies in the circle with center xyc and radius r

%Version 2017.1
%Copyright 2014-2017 Paolo Bardella
if exist('bsxfun','builtin')
    %Only for MATLAB>=2007a
    d=bsxfun(@minus,xy,xyc);   
else
    d=xy-repmat(xyc,size(xy,1),1);
end
 res=(d(:,1).^2+d(:,2).^2)<r*r;
