# sampledPlaneAverage

[![OpenFOAM version](https://img.shields.io/badge/OpenFOAM-7-brightgreen)](https://github.com/OpenFOAM/OpenFOAM-7)

## functions
- sampledSet/circle
- sampledSurface/sampledBoundedPlane
- sampledSurface/sampledPlaneAverage
- sampledSurface/sampledPlaneSpanwise

## usage
```
functions
{
    surfaceSampling1
    {
        type surfaces;
        libs ("libsampledPlaneAverage.so");
        enabled true;
        writeControl outputTime; // The same time control with case data write.
        interpolationScheme  cell; // suggest to use 'cell', which is suitable for various mesh size in the domain.
        surfaceFormat vtk;
        fields ( p T );
        surfaces
        (
            plane1
            {
                type sampledPlaneSpanwise; // choose class 'spanwiseAverage' to get spanwise averaging.
                interpolate    true; // Must be true, or it just cut the plane do not do spanwise average.
                planeType pointAndNormal; // method used to specific the cuttingPlane
                pointAndNormalDict
                {
                    basePoint (0 0 0.1); // This is the origin point for cut plane and for spanwise average
                    normalVector (1 0 0); // This is the normal vector for cut plane, it must be one of the axis.
                }
                nPoints    72; // The sample points, how many data in 360 degree is used to average.
                spanwiseVector (0 0 1); // The symmetry axis, it is the rotation vector for spanwise average
                                        // it must be perpendicular to cut plan normalVector
            }
        );
    }
    surfaceSampling2
    {
        type surfaces;
        libs ("libsampledPlaneAverage.so");
        enabled true;
        writeControl outputTime; // The same time control with case data write.
        interpolationScheme  cell; // suggest to use 'cell', which is suitable for various mesh size in the domain.
        surfaceFormat vtk;
        fields ( p T );
        surfaces
        (
            plane2
            {
                type    sampledPlaneAverage;
                planeType   pointAndNormal;
                pointAndNormalDict
                {
                    normal (1 0 0);
                    point  (0 0 0);
                }
                interpolate	true;		//must be true
                distance	1E-4;		//interval
                lineOfSight	false;		// false is average; true is sum
            }
        );
    }    
}
```
