function ism = ismatrix ( u ) ;

ism = (ndims(u) == 2) & (min(size(u)) ~= 1);
