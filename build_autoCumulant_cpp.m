% Please install "Microsoft Visual Studio 2012 Redistributable Version 4"
% for the compiled mex files to work correctly
%
% Author: Simon Christoph Stein
% E-Mail: scstein@phys.uni-goettingen.de
% Date: 2017

%%
% Build the cumulant file with openMP support (parallel processing)
% Using /MT statically links the visual studio runtime, but compiling with
% openMP creates a non-statically linkable runtime dependency. So with
% openMP the Microsoft Visual Studio 2012 Redistributable Version 4 is
% needed.

mex -v autoCumulant2D_cpp.cpp COMPFLAGS="$COMPFLAGS /openmp"
mex -v autoCumulant3D_cpp.cpp COMPFLAGS="$COMPFLAGS /openmp"

%% Without openMP (no parallel processing)
% mex -v autoCumulant2D_cpp.cpp COMPFLAGS="$COMPFLAGS /MT"
% mex -v autoCumulant3D_cpp.cpp COMPFLAGS="$COMPFLAGS /MT"
