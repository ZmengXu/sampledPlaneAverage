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

#include "sampledPlaneAverage.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPlaneAverage, 0);
    //addNamedToRunTimeSelectionTable(sampledSurface, sampledPlaneAverage, word, spatialAverage);// the last argument --spatialAverage,we use it in controlDict
    addToRunTimeSelectionTable(sampledSurface, sampledPlaneAverage, word);// we use sampledPlaneAverage in controlDict
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*Foam::sampledPlaneAverage::sampledPlaneAverage
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const scalar& distance,
    const keyType& zoneKey,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    cuttingPlane(planeDesc),
    distance_(distance),
	lineOfSight_(false),
    zoneKey_(zoneKey),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug && zoneKey_.size() && mesh.cellZones().findIndex(zoneKey_) < 0)
    {
        Info<< "cellZone " << zoneKey_
            << " not found - using entire mesh" << endl;
    }

    const vector& normal = planeDesc.normal()/mag(planeDesc.normal());

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
}*/


Foam::sampledPlaneAverage::sampledPlaneAverage
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    //cuttingPlane(plane(dict.lookup("basePoint"), dict.lookup("normalVector"))),
    cuttingPlane(plane(dict)),
    distance_(readScalar(dict.lookup("distance"))),
	lineOfSight_(dict.lookupOrDefault<bool>("lineOfSight", false)),
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

    const vector& normal = this->normal()/mag(this->normal());

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

Foam::sampledPlaneAverage::~sampledPlaneAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledPlaneAverage::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledPlaneAverage::expire()
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


bool Foam::sampledPlaneAverage::update()
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


Foam::tmp<Foam::scalarField> Foam::sampledPlaneAverage::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlaneAverage::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlaneAverage::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlaneAverage::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlaneAverage::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledPlaneAverage::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlaneAverage::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlaneAverage::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlaneAverage::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlaneAverage::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledPlaneAverage::print(Ostream& os) const
{
    os  << "sampledPlaneAverage: " << name() << " :"
        << "  base:" << refPoint()
        << "  normal:" << normal()
        << "  triangulate:" << triangulate_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
