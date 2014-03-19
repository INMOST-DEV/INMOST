/*
 * point.h
 *
 *  Created on: 13.09.2011
 *      Author: chernyshenko
 */

#ifndef POINT_H_
#define POINT_H_
/*********************************
 Description of points

 Information
 DIM         - space dimension
 COORD_TYPE  - type of point coordinats,
               must be comparable for
			   sort to work
 SCALE_TYPE  - type of index for array
               must be integrable

 Functions
 swap     - swap two points
 sort     - implimentation of quicksort
            perdimensional sort
 gen      - generate random set
 destroy  - release memory
 print    - write data to console
 step     - makes one step by updating
			velocities and coordinates
 draw     - draws set of points

 Dependency: point.h
   Standard: stdlib.h, time.h, math.h
   Specific: glut.h
**********************************/

#ifndef _POINT_H
#define _POINT_H

#define DIM          3
#define COORD_TYPE   float
#define SCALE_TYPE   unsigned int

struct point
{
	COORD_TYPE coord[DIM];
	COORD_TYPE radius[DIM];
	int polynum;
};

#endif


#endif /* POINT_H_ */
