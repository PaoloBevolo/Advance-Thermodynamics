Last updated: 10/03/2016

General
-------
The library uses two externally developped projects:
* the combination of shapes is based on "Polygon Clipper" by S. Holz
  http://www.mathworks.com/matlabcentral/fileexchange/8818-polygon-clipper

* the program TRIANGLE, A Two-Dimensional Quality Mesh Generator and 
  Delaunay Triangulator, by Jonathan Richard Shewchuk
  https://www.cs.cmu.edu/~quake/triangle.html
  modified to support the Microsoft Visual Studio Compiler 2015 and wrapped 
  in a mex file by Paolo Bardella.

How to recompile the binary files
---------------------------------

The library is distributed with compiled mex files which should work on 32 and 64 bit Windows and *nix systems and has been tested with MATLAB version 2015a and 2015b.
 
To recompile polygon-clipper:

cd('PolygonClipper');
mex gpc.c gpc_mexfile.c -output PolygonClip
movefile(['PolygonClip.' mexext()],fullfile('..','..','private',['PolygonClip.' mexext()]),'f')


To recompile Triangle (MEX version):

cd('TriangleMex');
mex trianglemex.c triangle.c -DTRILIBRARY -DWIN32 -DTRIANGLEMEX_EXPORTS;
movefile(['TriangleMEX.' mexext()],fullfile('..','..','private',['TriangleMEX.' mexext()]),'f')


