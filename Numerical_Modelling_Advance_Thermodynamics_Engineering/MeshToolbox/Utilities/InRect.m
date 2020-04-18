function res=InRect(xy,xyc,wh)
%Return if a point (x,y) lies in the rectangle with center xyc and size wh
% InRect([1,1],[0,0],[5,3]);
%Version 2017.1
%Copyright 2014-2017 Paolo Bardella

if exist('bsxfun','builtin')
    %Only for MATLAB>=2007a
    res=all(bsxfun(@lt,abs(bsxfun(@minus,xy,xyc)),wh/2)')';
else
    res=all((abs(xy-repmat(xyc,size(xy,1),1))-repmat(wh/2,size(xy,1),1))'<=0)';   
end