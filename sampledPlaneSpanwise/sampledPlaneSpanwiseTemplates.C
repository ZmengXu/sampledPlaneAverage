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

#include "sampledPlaneSpanwise.H"
#include "circleSet.H"
#include "meshSearch.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlaneSpanwise::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    return tmp<Field<Type>>(new Field<Type>(vField, meshCells()));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlaneSpanwise::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
	
    // One value per point.
    // Initialize with Zero to handle missed/degenerate faces

    tmp<Field<Type>> tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues.ref();

    boolList pointDone(points().size(), false);

    const faceList& fcs = faces();
    
    // Mesh search engine, no triangulation of faces.
    meshSearch searchEngine(mesh());//, polyMesh::FACE_PLANES
    
	// base Point and vector for cutPlane
	const point&  basePoint  = this->refPoint();
	const vector& baseVector = spanwiseVector_;
	const scalar dTheta = 360/nPoints_;//2*constant::pi/nPoints_;
	const vector normalVector = baseVector/mag(baseVector);

	if(debug)
	{
		Info << "basePoint baseVector dTheta normalVector and plane vector" << basePoint <<baseVector<< dTheta <<endl<< normalVector<<endl<< this->normal() <<endl;
	}

    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];
        for (const label pointI : f)
        {       
            if (!pointDone[pointI])
            {
                //- Creation of the a circle line that start from points()[pointI] with nPoints_ to form a circle.
                const point& currentPoint = points()[pointI];

				if(debug)
				{
					Info << " 01 " << " pointI " << pointI << " currentPoint " << currentPoint << endl;
				}
					
				vector circleVector  = currentPoint - basePoint;
				point  circleOrigin  = basePoint + ( circleVector & normalVector ) * normalVector;
				vector circleAxis    = currentPoint - circleOrigin;
				//vector circleAxis = (currentPoint - circleOrigin)/mag(currentPoint - circleOrigin);


				if(debug)
				{
					Info << " 02 " << " circleOrigin " << circleOrigin << " circleVector " << circleVector << " circleAxis " << circleAxis << endl;
				}

				if( mag(circleVector) < SMALL || mag(circleAxis) < SMALL)
				{
					values[pointI] = interpolator.interpolate
					(
						points()[pointI],
						meshCells()[faceI]
					);
				}
				else
				{
					tmp<Field<Type> > tlinevalues(new Field<Type>(nPoints_));
					Field<Type>& linevalues = tlinevalues.ref();

					// Extract the circle in the homogeneous direction 
					circleSet	line("circleLine", mesh(), searchEngine, axis_, circleOrigin, normalVector, currentPoint, dTheta);

					if(debug)
					{
						Info << " 03 " << " line.size() " << line.size() << endl;
					}

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
					
					if(debug)
					{
						Info << " 04 " << " tlinevalues " << tlinevalues.ref() << endl;
					}
					
					//- Compute the average of the line
					values[pointI] = Foam::average(tlinevalues);
					
				}
			    
				if(debug)
				{
					Info << "Average value : " << values[pointI] << endl;
				}
				
				//- Register that the point was analyzed  
                pointDone[pointI] = true;
            }
        }
    }

    return tvalues;

}


// ************************************************************************* //

