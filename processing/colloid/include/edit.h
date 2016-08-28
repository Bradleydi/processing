/*!
 * Dynamical properties.
 */

#ifndef EDIT_H
#define EDIT_H

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
//void gr2d( colloid_base&,
//		   const colloid_base&,
//		   bool =true);

void pt_qchi4(colloid_base&, double , int =2);

void ptid_qchi4(colloid_base&, double , int =2);
void block_ptid_qchi4(colloid_base&, double , int =2);

#endif /* EDIT_H */
