function RebuildMex()
%rebuild the MEX files used in the library
%Version 2017.1
%Copyright 2014-2017 Paolo Bardella

clear('mex');                                                          %#ok
CurrentFolder=pwd;
FunctionFolder=fileparts(mfilename('fullpath'));
%first, compile PolygonClipper 
try
    disp('Building PolygonClipper...');
    cd(fullfile(FunctionFolder,'PolygonClipper'));
    mex gpc.c gpc_mexfile.c -output PolygonClip
    movefile(['PolygonClip.' mexext()],fullfile('..','..','private',['PolygonClip.' mexext()]),'f')
    disp('Done.');
catch Exc
    warning('An error occurred building PolygonClipper');
    disp(Exc);
end
%% then compile TriangleMEX
try
    disp('Building TriangleMEX...');
    cd(fullfile(FunctionFolder,'TriangleMEX'));
    mex TriangleMEX.c triangle.c -DTRILIBRARY -DWIN32 -DTRIANGLEMEX_EXPORTS;
    movefile(['TriangleMEX.' mexext()],fullfile('..','..','private',['TriangleMEX.' mexext()]),'f')
    disp('Done.');
catch Exc
    warning('An error occurred building TriangleMEX');
    disp(Exc);
end

cd(CurrentFolder);