// version 0.9


/*!
 * Static structural properties.
 */

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "colloid_base.h"

/*! 2D pair correlation function (radial distribution function)
 *  @param gr
 *  	pair correlation function result
 *  @param p
 *  	position data (tracked or not)
 *  @param tracked
 *  	A tracked data or not:
 *  	true:	tracked, defalt case
 *  	false:	not tracked
 *
 *  I'm really lazy that I want to avoid calculating the arc area interseted by
 *  the rectangular field of view box. So I restrict the particles accounted as
 *  the ones that all particles within the maximum radius to calculated lies 
 *  in the field of view.
 *
 *  Actually, the structural properties can be extracted from the untracked
 *  data, but my code is not good now, since untracked data could not be used
 *  now. I should try to solve this program recently, however, before that I
 *  should make sure the program now is correct.
 */
void gr2d(const colloid_base&, bool =true);


/*! PURPOSE: 
 *  	bond orientational order parameter
 *  	psi_4: square lattice
 *  	psi_6: triangular lattice
 *
 *  	\psi_n=\sum_{(i,j)}\exp(in\theta_{ij}),
 *  	
 *  	where
 *  		the sum is over all $i$'s nearest neighbor $j$,
 *  		and $theta_{ij}$ is the angle of the bond $ij$.
 *
 *  @param p
 *  	(type) colloid_base&
 *  	(input) data to calculate the orientational order parameter.
 *  @param n
 *  	(type) int
 *  	(input) nearest neighbor number.
 *  @param Re_psi
 *  	(type) double *
 *  	(output) real part of psi.
 *  @param Im_psi
 *  	(type) double *
 *  	(output) imaginary part of psi.
 *  @param tracked
 *  	(type) bool
 *  	(input) whether input data is tracked.
 *  	(default) true, tracked.
 *
 *	WARNING: 
 *  	Re_psi and Im_psi should be allocated memory by
 *  		(double *)malloc((p.get_size()[1])*sizeof(double));
 * 		before calling the function.
 *
 * 	This function calls dct or dct_RNG defined in dct.h, which is used for
 * 	the triangulation and determining the nearest neighbors for each particle.
 */
void psi_n(colloid_base &, int, double *, double *, bool =true);


/*! PURPOSE: 
 *  	draw bond orientational order parameter
 *  	psi_4: square lattice
 *  	psi_6: triangular lattice
 *
 *  	\psi_n=\sum_{(i,j)}\exp(in\theta_{ij}),
 *  	
 *  	where
 *  		the sum is over all $i$'s nearest neighbor $j$,
 *  		and $theta_{ij}$ is the angle of the bond $ij$.
 *
 *  @param p
 *  	(type) colloid_base&
 *  	(input) data to calculate the orientational order parameter.
 *  @param n
 *  	(type) int
 *  	(input) nearest neighbor number.
 *  @param filename
 *  	(type) const char *
 *  	(input) fig filename.
 *  @param tracked
 *  	(type) bool
 *  	(input) whether input data is tracked.
 *  	(default) true, tracked.
 *
 *	WARNING: 
 *  	Re_psi and Im_psi should be allocated memory by
 *  		(double *)malloc((p.get_size()[1])*sizeof(double));
 * 		before calling the function.
 *
 * 		Since the function outputs figures in fig format, which can be
 * 		handled by Xfig, so please be happy to use 
 * 			fig2dev -L eps -g white file.fig file.eps
 * 		to convert the fig file to be an eps figure, with white background.
 *
 * 		Please be happy, since Xfig is quite nice.
 *
 * 	This function calls dct or dct_RNG defined in dct.h, which is used for
 * 	the triangulation and determining the nearest neighbors for each particle.
 */
void draw_psi_n(colloid_base &, int, const char *, bool =true);


void draw_psi_n_with_bond(colloid_base &, int ,
        const char *, int =100, bool =true);

void psi_n_statistics(colloid_base &, int , 
		double =0.5, bool =true);

void triangle(const char *file, double threshold, bool tracked);
void triangle_config(const char *file, double threshold, bool tracked);
void triangle_cluster(const char *file, double threshold, bool tracked);
void triangle_cluster_plot(const char *file, int frame, bool tracked);

void plot_nucleus(const char *file, const char *str_frameindex, 
				  int framerate, bool tracked);

void s2(const char *file, const int dim, bool tracked);

void sk(const char *file, bool Sk2D, bool tracked);
#endif /* STRUCTURE_H */		 
