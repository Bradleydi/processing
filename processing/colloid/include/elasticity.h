#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "colloid_base.h"

/*! calculate the partial derivatives of displacement from tracked trajectory
 * @param cg
 * 		store the coase-grained partial derivatives of displacment
 * 		{center_x, center_y, ..., u_xx, u_xy, ..., u_yx, u_yy, ..., ..., Dx^2, Dy^2, ..., D^2, t}
 * 		(dim*dim+2*dim+2, total_box_in_all_frame, data)
 * @param ptid
 * 		tracked position data
 * @param b
 * 		coarse grain constant, L_b=L/b
 * @param Lambda
 * 		if the distance of two particles is smaller than Lambda, then 
 * 		the two particles are neighbors of each other.
 * @param Remove_Drift
 * 		true: remove the drift; false: don't remove drift
 */
void coarse_grain(colloid_base&, 
		          const colloid_base&, 
				  const int&, 
				  const float &, 
				  bool =true);

void elasticity(const colloid_base&,
		        const int &,
				const float &,
				const int &,
				const char * ,
				const float & =0.1,
				bool =true);

void el_modulus(const char*file, int b, double a0, int dim);


void strain( const char *file, 
				  const float Lambda,
				  bool Remove_Drift);


void strain_plot(const char *file, const char *str_frameindex, 
				const char *which_strain);

void plot_strain(const char *file, const char *str_frameindex, 
				const char *which_strain, int framerate);

#endif /* ELASTICITY_H */
