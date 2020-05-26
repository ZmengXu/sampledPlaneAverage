/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "uniformSet.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlaneAverage::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    return tmp<Field<Type>>(new Field<Type>(vField, meshCells()));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlaneAverage::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per point.
    // Initialize with Zero to handle missed/degenerate faces

    tmp<Field<Type>> tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues.ref();

	// Get the boundary points
	DynamicList<point> boundaryPoints;
	boundaryPoints.clear();
    const pointField& meshPoints = mesh().points();
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (!pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                const face& f = mesh().faces()[facei++];
                forAll(f, fp)
                {
                    boundaryPoints.append( meshPoints[f[fp]] );
                }
            }
        }
    }

	const pointField& cuttPoints = points();
    // Determine if point is on boundary. Points on boundaries will
	// not Extract the line in the homogeneous direction
    // Coupled boundaries are handled explicitly so not marked here.
    PackedBoolList isBoundaryPoint(cuttPoints.size(), false);

	forAll(cuttPoints, pointI)
	{
		forAll(boundaryPoints, boundaryPointI)
		{
			if( mag(cuttPoints[pointI]-boundaryPoints[boundaryPointI]) < SMALL )
			{
				isBoundaryPoint[pointI] = true;
			}
		}
	}

    boolList pointDone(cuttPoints.size(), false);
    const faceList& fcs = faces();
    
    // Mesh search engine, no triangulation of faces.
    meshSearch searchEngine(mesh());//, polyMesh::FACE_PLANES
    
	// base Point and vector for cutPlane
	const point&  basePoint  = this->refPoint();
	const vector normalVector = this->normal()/mag(this->normal());

	const scalar dInterval = distance_;

	const vectorField& originalCenter(mesh().cellCentres());
	const scalarField projectedDistan(((originalCenter - basePoint) & normalVector));

	const label startIndex = min(projectedDistan)/dInterval;
	const label endIndex   = max(projectedDistan)/dInterval;
	const label nPoints    = endIndex - startIndex;

	if( debug >= 1 )
	{
		Info << "basePoint dInterval and plane vector" << basePoint << "	" << dInterval<< "	" << normalVector <<endl;
		Info << " 11 "  << "startIndex " << startIndex << "	endIndex " << endIndex << "	nPoints " << nPoints << endl;
	}

    forAll(fcs, faceI)
    {
		if( debug >= 2 )
		{
			Info << " 20 " << " faceI " << faceI << " totalFace " << fcs.size() << endl;
		}
        const face& f = fcs[faceI];
        for (const label pointI : f)
        {
            if ( !pointDone[pointI] )
            {
                //- Creation of the a straight line that start from points()[pointI] with nPoints_ to form a line.
                const point& currentPoint = points()[pointI];

				if( debug >= 3 )
				{
					Info << " 31 " << " pointI " << pointI << " currentPoint " << currentPoint << endl;
				}

				point startPoint = currentPoint + startIndex*dInterval*normalVector;
				point endPoint   = currentPoint +   endIndex*dInterval*normalVector;

				if( debug >= 3 )
				{
					Info << " 32 " << " startPoint " << startPoint << " endPoint " << endPoint << " nPoints " << nPoints << endl;
				}

				if(!isBoundaryPoint[pointI])
				{
					// Extract the line in the homogeneous direction 	
					uniformSet line("homogeneousLine", mesh(), searchEngine, axis_, startPoint, endPoint, nPoints);

					if( line.size() == 0 )
					{
						if( debug >= 2 )
						{
							Info << " 21 " << " uniformSet line has a size of " << line.size() << '\n'
									<< " The current point line can not found " << currentPoint << endl;
						}
						values[pointI] *= 0;
					}
					else
					{
						if( debug >= 3 )
						{
							Info << " 33 " << " line.size() " << line.size() << endl << " line " << line << endl;
						}

						Field<Type> linevalues(line.size());

						//- Interpolate the values of the fields in each point of the line
						forAll(line, lineI)
						{
							linevalues[lineI] = interpolator.interpolate
							(
								line[lineI],
								line.cells()[lineI],
								line.faces()[lineI]
							);
						}
				
						if( debug >= 3 )
						{
							Info << " 34 " << " linevalues " << linevalues << endl;
						}
						if(lineOfSight_)
						{
							values[pointI] = Foam::sum(linevalues);
						}
						else
						{
							//- Compute the average of the line
							values[pointI] = Foam::average(linevalues);
						}

					}
				}
				else
				{
					if( debug >= 2 )
					{
						Info << " 22 " << " pointI " << pointI << " is on the boundary face " << currentPoint << endl;
					}
					values[pointI] *= 0;
				}

				if( debug >= 3 )
				{
					Info << " 35 "  << "Average value : " << values[pointI] << endl;
				}
			
				//- Register that the point was analyzed  
                pointDone[pointI] = true;
            }
        }
    }

    return tvalues; 
}


// ************************************************************************* //
