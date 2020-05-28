/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
namespace sampledSurfaces
{
    defineTypeNameAndDebug(sampledPlaneSpanwise, 0);
    addToRunTimeSelectionTable(sampledSurface, sampledPlaneSpanwise, word);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*Foam::sampledSurfaces::sampledPlaneSpanwise::sampledPlaneSpanwise
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
    Info<<"I am in sampledPlaneSpanwise::sampledPlaneSpanwise1 =========sampledPlaneSpanwise======="<<endl;

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
}*/


Foam::sampledSurfaces::sampledPlaneSpanwise::sampledPlaneSpanwise
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    //cuttingPlane(Foam::plane(dict)),
    cuttingPlane(plane(dict)),
    nPoints_(dict.lookupOrDefault<label>("nPoints", 12)),
	spanwiseVector_(vector(dict.lookup("spanwiseVector"))),
    zoneKey_(keyType::null),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    Info<<"I am in sampledPlaneSpanwise::sampledPlaneSpanwise2 =========sampledPlaneSpanwise======="<<endl;

    // Make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found("coordinateSystem"))
    {
        coordinateSystem cs(mesh, dict.subDict("coordinateSystem"));

        point  base = cs.globalPosition(planeDesc().refPoint());
        vector norm = cs.globalVector(planeDesc().normal());

        // Assign the plane description
        static_cast<Foam::plane&>(*this) = Foam::plane(base, norm);
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
    Pout<< "axis_==="<<axis_<<endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledPlaneSpanwise::~sampledPlaneSpanwise()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::sampledPlaneSpanwise::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledSurfaces::sampledPlaneSpanwise::expire()
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


bool Foam::sampledSurfaces::sampledPlaneSpanwise::update()
{
    Info<<"I am in sampledPlaneSpanwise::update =========sampledPlaneSpanwise======="<<endl;

    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    labelList selectedCells = mesh().cellZones().findMatching(zoneKey_).used();

    bool fullMesh = returnReduce(selectedCells.empty(), andOp<bool>());

    if (fullMesh)
    {
        const label len = mesh().nCells();

        selectedCells.setSize(len);

        label count = 0;
        for (label celli=0; celli < len; ++celli)
        {
            selectedCells[count++] = celli;
        }
    }
    else
    {
        label count = 0;
        for (const label celli : selectedCells)
        {
            selectedCells[count++] = celli;
        }

        selectedCells.setSize(count);
    }

    reCut(mesh(), triangulate_, selectedCells);

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledSurfaces::sampledPlaneSpanwise::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledSurfaces::sampledPlaneSpanwise::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurfaces::sampledPlaneSpanwise::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurfaces::sampledPlaneSpanwise::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledSurfaces::sampledPlaneSpanwise::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledSurfaces::sampledPlaneSpanwise::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledSurfaces::sampledPlaneSpanwise::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurfaces::sampledPlaneSpanwise::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurfaces::sampledPlaneSpanwise::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledSurfaces::sampledPlaneSpanwise::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledSurfaces::sampledPlaneSpanwise::print(Ostream& os) const
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
