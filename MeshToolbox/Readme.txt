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
Detailed instructions are provided in the Readme.txt file in the 
OriginalDependencies folder.

A script, named RebuildMex.m, can be used to rebuild all the dependencies and 
to copy the generated files to the correct folder
