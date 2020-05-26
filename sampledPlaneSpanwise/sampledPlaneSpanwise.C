/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sampledPlaneSpanwise.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPlaneSpanwise, 0);
    addNamedToRunTimeSelectionTable(sampledSurface, sampledPlaneSpanwise, word, spanwiseAverage);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPlaneSpanwise::sampledPlaneSpanwise
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
	const vector& spanwiseVector,
    const keyType& zoneKey,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    cuttingPlane(planeDesc),
    nPoints_(12),
	spanwiseVector_(spanwiseVector),
    zoneKey_(zoneKey),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
  
    const vector& normal = spanwiseVector_;
    const vector& planNorm = planeDesc.normal();
    if ( abs(normal&planNorm) > SMALL)
    {
        FatalErrorIn
	    (
	         "Foam::sampledAveragePlane::sampledAveragePlane"
	    )   << "The spanwiseVector normal is not perpendicular to the plane normal"
	        << exit(FatalError);
    }
    if (abs(abs(normal.x())+abs(normal.y())+abs(normal.z())-1.0) > SMALL)
    {
        FatalErrorIn
	    (
	         "Foam::sampledAveragePlane::sampledAveragePlane"
	    )   << "The spanwiseVector normal is not one of the coordinate axis"
	        << exit(FatalError);
    }
    
    if(abs(abs(normal.x())-1.0) < SMALL)
    {
        axis_ = "x";
    }
    else if(abs(abs(normal.y())-1.0) < SMALL)
    {
        axis_ = "y";
    }
    else if(abs(abs(normal.z())-1.0) < SMALL)
    {
        axis_ = "z";
    }
}


Foam::sampledPlaneSpanwise::sampledPlaneSpanwise
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    cuttingPlane(plane(dict.lookup("basePoint"), dict.lookup("normalVector"))),
    nPoints_(dict.lookupOrDefault<label>("nPoints", 12)),
	spanwiseVector_(vector(dict.lookup("spanwiseVector"))),    
    zoneKey_(keyType::null),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    // make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(mesh, dict.subDict("coordinateSystem"));

        point  base = cs.globalPosition(planeDesc().refPoint());
        vector norm = cs.globalVector(planeDesc().normal());

        // assign the plane description
        static_cast<plane&>(*this) = plane(base, norm);
    }

    dict.readIfPresent("zone", zoneKey_);

    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }
   
    const vector& normal = spanwiseVector_;
    const vector& planNorm = this->normal();
	if ( mag(normal&planNorm) > SMALL)
    {
        FatalErrorIn
	    (
	         "Foam::sampledAveragePlane::sampledAveragePlane"
	    )   << "The spanwiseVector normal is not perpendicular to the plane normal"
	        << exit(FatalError);
    }
    if (abs(abs(normal.x())+abs(normal.y())+abs(normal.z())-1.0) > SMALL)
    {
        FatalErrorIn
	    (
	         "Foam::sampledAveragePlane::sampledAveragePlane"
	    )   << "The spanwiseVector normal is not one of the coordinate axis"
	        << exit(FatalError);
    }
    
    if(abs(abs(normal.x())-1.0) < SMALL)
    {
        axis_ = "x";
    }
    else if(abs(abs(normal.y())-1.0) < SMALL)
    {
        axis_ = "y";
    }
    else if(abs(abs(normal.z())-1.0) < SMALL)
    {
        axis_ = "z";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledPlaneSpanwise::~sampledPlaneSpanwise()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledPlaneSpanwise::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledPlaneSpanwise::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledPlaneSpanwise::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    labelList selectedCells = mesh().cellZones().findMatching(zoneKey_).used();

    if (selectedCells.empty())
    {
        reCut(mesh(), triangulate_);
    }
    else
    {
        reCut(mesh(), triangulate_, selectedCells);
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledPlaneSpanwise::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlaneSpanwise::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlaneSpanwise::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlaneSpanwise::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlaneSpanwise::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledPlaneSpanwise::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlaneSpanwise::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlaneSpanwise::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlaneSpanwise::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlaneSpanwise::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledPlaneSpanwise::print(Ostream& os) const
{
    os  << "sampledPlaneSpanwise: " << name() << " :"
        << "  base:" << refPoint()
        << "  normal:" << normal()
        << "  nPoints:" << nPoints_        
        << "  triangulate:" << triangulate_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
template <class Type>
Foam::tmp<Foam::Tensor<Type>>
Foam::sampledPlaneSpanwise::rotationMatrix(vector basePoint, vector baseVector, scalar theta)
{
    if (debug)
    {
        Pout<< "sampledSpanwise::rotationMatrix : Create an matrix to ratate."
            << endl;
    }

	tmp<Tensor<Type>> tTransformationMatrix(new Tensor<Type>(0.0));    
    Tensor<Type>& TransformationMatrix = tTransformationMatrix();

	scalar a = basePoint.x();
	scalar b = basePoint.y();
	scalar c = basePoint.z();

	vector normalVector = baseVector/mag(baseVector);
	scalar u = normalVector.x();
	scalar v = normalVector.y();
	scalar w = normalVector.z();

	scalar uu = u * u;
	scalar uv = u * v;
	scalar uw = u * w;
	scalar vv = v * v;
	scalar vw = v * w;
	scalar ww = w * w;
	scalar au = a * u;
	scalar av = a * v;
	scalar aw = a * w;
	scalar bu = b * u;
	scalar bv = b * v;
	scalar bw = b * w;
	scalar cu = c * u;
	scalar cv = c * v;
	scalar cw = c * w;

	scalar costheta = cos(theta);
	scalar sintheta = sin(theta);

	TransformationMatrix[0][0] = uu + (vv + ww) * costheta;
	TransformationMatrix[0][1] = uv * (1 - costheta) + w * sintheta;
	TransformationMatrix[0][2] = uw * (1 - costheta) - v * sintheta;
	TransformationMatrix[0][3] = 0;

	TransformationMatrix[1][0] = uv * (1 - costheta) - w * sintheta;
	TransformationMatrix[1][1] = vv + (uu + ww) * costheta;
	TransformationMatrix[1][2] = vw * (1 - costheta) + u * sintheta;
	TransformationMatrix[1][3] = 0;

	TransformationMatrix[2][0] = uw * (1 - costheta) + v * sintheta;
	TransformationMatrix[2][1] = vw * (1 - costheta) - u * sintheta;
	TransformationMatrix[2][2] = ww + (uu + vv) * costheta;
	TransformationMatrix[2][3] = 0;

	TransformationMatrix[3][0] = (a * (vv + ww) - u * (bv + cw)) * (1 - costheta) + (bw - cv) * sintheta;
	TransformationMatrix[3][1] = (b * (uu + ww) - v * (au + cw)) * (1 - costheta) + (cu - aw) * sintheta;
	TransformationMatrix[3][2] = (c * (uu + vv) - w * (au + bv)) * (1 - costheta) + (av - bu) * sintheta;
	TransformationMatrix[3][3] = 1;


    if (debug)
    {
        Pout<< "TransformationMatrix"	<< endl
			<< "xx xy xz" << TransformationMatrix.xx() << TransformationMatrix.xy() << TransformationMatrix.xz() << endl
			<< "yx yy yz" << TransformationMatrix.yx() << TransformationMatrix.yy() << TransformationMatrix.yz() << endl
			<< "zx zy zz" << TransformationMatrix.zx() << TransformationMatrix.zy() << TransformationMatrix.zz() << endl        
            << "Transformation from basic point and normal vector"	<< endl;
    }

	return tTransformationMatrix;

}
