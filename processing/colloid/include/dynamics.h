/*!
 * Dynamical properties.
 */

#ifndef DYNAMICS_H
#define DYNAMICS_H

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
void dxpdf(colloid_base&, const int & =1, const int & =2, const int & =100);

void msd(colloid_base&, const int & =2, bool =true);

void DWfactor(const char *, double =-1, const int & =2, bool =true);
//void qchi4(colloid_base&, float&);

void vibration_Gaussian(const char *file, bool Remove_drift);

#endif /* DYNAMICS_H */
