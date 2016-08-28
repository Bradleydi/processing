#ifndef DATA_PREPROCESS_H
#define DATA_PREPROCESS_H

#include "colloid_base.h"

/*! select trajectories longer than a given number
 * @param cb
 * 		to store the result
 * @param ptid
 * 		(postion, t, id) data
 * @param n
 * 		least trajector length to select
 */
void trajectory_filter(colloid_base &, const colloid_base &, const int &);

/*! select trajectories longer than a given number
 * @param outfile
 * 		to store the result
 * @param infile
 * 		file containing (postion, t, id) data
 * @param n
 * 		least trajector length to select
 */
void trajectory_filter(const char *, const char *, const int &);

/*! select trajectories longer than a given number, the output file is
 *  automatically created as
 *      infile(remove ".gdf")_n.gdf
 * @param infile
 * 		file containing (postion, t, id) data
 * @param n
 * 		least trajector length to select
 */
void trajectory_filter(const char *, const int &);


/*! select data in spatial region 
 *      ( xmin_new, xmax_new ) x ( ymin_new, ymax_new )
 *  and select trajectories longer than a given number
 * @param cb
 * 		to store the result
 * @param ptid
 * 		(postion, t, id) data
 * @param xmin_new
 *      the new xmin
 * @param xmax_new
 *      the new xmax
 * @param ymin_new
 *      the new ymin
 * @param ymax_new
 *      the new ymax
 * @param n
 * 		least trajector length to select
 */
void space_filter(colloid_base &, 
		          const colloid_base &,
				  const float & ,
				  const float & ,
				  const float & ,
				  const float & ,
				  const int & );


/*! remove the drift for pt data
 * @param pt
 * 		data to remove the drift, results will store in it.
 *
 */
void pt_remove_drift(colloid_base& );

void ptid2pt(colloid_base &, colloid_base &);

void ptid_drift(const char *, int =0);

void remove_drift(colloid_base &, int =0);


void frame_cut(const char *, const int, bool =false);

void showposition(const char *, const int =1, bool =false);

void space_cut(const char *, const char *);


void traj_length(const char *file);


void track(const char *file, const float maxdisp, const int memory);

bool is_pt(const char *file);
bool is_ptid(const char *file);

void goodtrack(const char *file, double threshold);

/* local linear regression of the drift */
void drift_llr(const char *file, int dim, int smooth);


void IR_SSD(const char *file, int dim, int window);
#endif /* DATA_PREPROCESS_H */
