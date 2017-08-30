## Fourier Super-resolution Optical Fluctuation Imaging (fSOFI)


MATLAB (& MEX/C++) implementation of Fourier Super-resolution Optical Fluctuation Imaging (fSOFI) as described in the publication:
	
 > Stein, S.C.; Huss, A.; Hähnel, D.; Gregor, I.; Enderlein, J.  
 >  [Fourier interpolation stochastic optical fluctuation imaging](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-23-12-16154),  
 >  Optics Express  23, 16154-16163, 2015.
 
 Stochastic Optical Fluctuation Imaging (SOFI) is a super-resolution fluorescence microscopy technique which allows to enhance the spatial resolution of an image by evaluating the temporal fluctuations of blinking fluorescent emitters. SOFI is not based on the identification and localization of single molecules such as in the widely used Photoactivation Localization Microsopy (PALM) or Stochastic Optical Reconstruction Microscopy (STORM), but computes a superresolved image via temporal cumulants from a recorded movie. A technical challenge hereby is that, when directly applying the SOFI algorithm to a movie of raw images, the pixel size of the final SOFI image is the same as that of the original images, which becomes problematic when the final SOFI resolution is much smaller than this value.

The Fourier SOFI (fSOFI) algorithm pre-processes the raw input movie to generate physically meaningful higher-resolution output images, yielding an improved image fidelity. While the function `FSOFI_Analysis.m` works on `2D+time` data, its counterpart `FSOFI_Analysis3D.m` processes `3D+time` data (e.g. acquired by recording multiple focal planes simultaneously using a prism).  
Note: In the latter case often only very few planes are available axially and it is strongly recommended to set the `mirrorMode='axial'` option, which prevents artifacts from non-periodic boundaries during Fourier interpolation.

It is recommended to use the high-performance MEX-functions for the computation of cumulants. To use the precompiled files you must install the `Visual Studio 2012 C++ Redistributable (Update 4)` available on the [Microsoft website](https://www.microsoft.com/de-de/download/details.aspx?id=30679). If a working Visual Studio compiler is configured in MATLAB, the C++ files can also be compiled using the included `build_\*.m` script. If, for some reason, MEX files are not an option, you can switch to equivalent (but slow) MATLAB implementations by uncommenting the respective lines in `FSOFI_Analysis.m` and `FSOFI_Analysis3D.m` (search for the occurrence of %*--- Core SOFI Computation ---*).