/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "uniformSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniformSet, 0);
    addToRunTimeSelectionTable(sampledSet, uniformSet, word);

    const scalar uniformSet::tol = 1e-3;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::uniformSet::nextSample
(
    const point& currentPt,
    const vector& offset,
    const scalar smallDist,
    point& samplePt,
    label& sampleI
) const
{
    bool pointFound = false;

    const vector normOffset = offset/mag(offset);

    samplePt += offset;
    sampleI++;

    for (; sampleI < nPoints_; sampleI++)
    {
        scalar s = (samplePt - currentPt) & normOffset;

        if (s > -smallDist)
        {
            // samplePt is close to or beyond currentPt -> use it
            pointFound = true;

            break;
        }
        samplePt += offset;
    }

    return pointFound;
}


void Foam::uniformSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    label sampleSegmentI = 0, sampleI = 0;
    scalar sampleT = 0;

    while (sampleI < nPoints_)
    {
        const point pt = (1 - sampleT)*start_ + sampleT*end_;

        const label sampleCellI = searchEngine().findCell(pt);

        if (sampleCellI == -1)
        {
            if (++ sampleI < nPoints_)
            {
                sampleT = scalar(sampleI)/(nPoints_ - 1);
            }
        }
        else
        {
            passiveParticle sampleParticle(mesh(), pt, sampleCellI);

            do
            {
                samplingPts.append(sampleParticle.position());
                samplingCells.append(sampleParticle.cell());
                samplingFaces.append(-1);
                samplingSegments.append(sampleSegmentI);
                samplingCurveDist.append(sampleT*mag(end_ - start_));

                if (++ sampleI < nPoints_)
                {
                    sampleT = scalar(sampleI)/(nPoints_ - 1);
                    sampleParticle.track((end_ - start_)/(nPoints_ - 1), 0);
                }
            }
            while (sampleI < nPoints_ && !sampleParticle.onBoundaryFace());

            ++ sampleSegmentI;
        }
    }
}


void Foam::uniformSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformSet::uniformSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const point& start,
    const point& end,
    const label nPoints
)
:
    sampledSet(name, mesh, searchEngine, axis),
    start_(start),
    end_(end),
    nPoints_(nPoints)
{
    genSamples();

    if (debug)
    {
        write(Pout);
    }
}


Foam::uniformSet::uniformSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    start_(dict.lookup("start")),
    end_(dict.lookup("end")),
    nPoints_(readLabel(dict.lookup("nPoints")))
{
    genSamples();

    if (debug)
    {
        write(Pout);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformSet::~uniformSet()
{}


// ************************************************************************* //
