#ifndef DCT_H
#define DCT_H

#include "DCT.h"

/*! divide-and-conquer triangulation
 * @param np
 * 		(input) number of points.
 * @param points
 * 		(input) points.
 * @param nedge
 * 		(output) number of edges.
 * @param edges
 * 		(output) edges, endpoints are indexed by index in points.
 * 		It should be allocated enough memory at the beginning.
 *		The maximum nedge=3*N,  
 *		one inner Delaunay triangle shares 3 edges with others, so
 *		3F=2E-b, where b is the edges on the convex hull.
 *		By Euler's characteristic V-E+F=1, for plane
 *		we have V-E+(2E-b)/3=1, so E=3(V-1)-b <=3V.
 */
void dct(int &, float *, int &, int *);


/*! divide-and-conquer relative neighborhood graph of spatial points,
 *  calculated by Delaunay triangulation.
 *  This function may be useful for determining neighbors in square lattice.
 * @param np
 * 	   (input) number of points.
 * @param points
 * 		(input) points.
 * @param nedge
 *	    (output) number of edges.
 * @param edges
 *      (output) edges, endpoints are indexed by index in points.
 *      It should be allocated enough memory at the beginning.
 *      The maximum nedge=3*N,  
 *      one inner Delaunay triangle shares 3 edges with others, so
 *      3F=2E-b, where b is the edges on the convex hull.
 *      By Euler's characteristic V-E+F=1, for plane
 *      we have V-E+(2E-b)/3=1, so E=3(V-1)-b <=3V.
 */
void dct_RNG(int &, float *, int &, int *);


/*! export resulting triangulaton or RNG as a fig file, 
 *  which can be manipulated by Xfig.
 * @param np
 * 	   (input) number of points.
 * @param points
 * 		(input) points.
 * @param nedge
 *	    (input) number of edges.
 * @param edges
 *      (input) edges, endpoints are indexed by index in points.
 *      edges should be calculated before calling this function.
 * @param filename
 * 		fig filename.
 * @param defect
 * 		defect type, 57 (triangular lattice) or 35 (square lattice),
 * 		default as 57.
 */
void dct_draw(int &, float *, int &, int *, const char *, unsigned char =57);

/*! return triangles from Delaunay triangulation
 * poiner tri shhould be allocated memory before calling this function
 * one inner Delaunay triangle shares 3 edges with others, so
 * 3F=2E-b, where b is the edges on the convex hull.
 * By Euler's characteristic V-E+F=1, for plane
 * we have V-(3F+b)/2+F=1, so F=2V-b <= 2V
 * so we can allocate tri by 6*V*sizeof(int)
 */
void dct_tri(int &np, float *points, int &ntri, int *tri);
#endif /* DCT_H */
