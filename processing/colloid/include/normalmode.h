#ifndef NORMAL_MODE_H
#define NORMAL_MODE_H


#include "colloid_base.h"
#include "lapack.h"
#include "statistics.h"

/* binargy io */

/*! write and read .ev file */
void writeev(int &, double *, double *, const char *);
int readev(int &, double **, double **, const char *);
/*! read only a mode from .ev file */
int readmode(int &, int &, double &, double **, const char *);
/*! read only eigenvalues from .ev file */
void readE(int &, double **,const char *);

/*! write and read .dcm file 
 * dynamic covariance matrix
 */
void writedcm(int &, double *, const char *);
void readdcm(int &, double **, const char *);

void writedcm3D(int &, double *, const char *);
void readdcm3D(int &, double **, const char *);
void dcm3D(const char *file);
void DOS3D(const char *file);


/*! ==================
 *  FORMAT of .mp FILE
 *  ==================
 *  header0 header1 mp
 *
 *  header0: int, 'm'*256+'p'
 *  header1: int, numbers of points * 2
 *  mp[header1] : float, xy positions
 *  	mp[2*i]   : x position of point i
 *  	mp[2*i+1] : y position of point i
 */
void writemp(const int &, float *, const char *);
int readmp(int &, float **, const char *);

/*! normal mode core 
 *  @param ptid
 *  	(type)    colloid_base &
 *  	(input)   [position, ...,  time, id] data
 *  @param file
 *  	(type)    const char *
 *  	(input)   filename prefix (can contain subfix) to store eigenvalues
 *  	          and eigenvectors of covariance matrix
 *  @param Remove_Drift
 *  	(type)    bool
 *  	(input)   remove drift or not
 *  	(default) remove drift
 *
 *  (input)ptid should contain same number of particles for each frame. 
 *  .ev file with prefix (input)file is generated, to store eigenvalues and
 *  eigenvectors of covariance matrix. The format is
 *  ==================
 *  FORMAT of .ev FILE
 *  ==================
 *  header0 header1 E G
 *
 *  header0: int, 'e'*256+'v'
 *  header1: int, numbers of eigenvalues
 *  E[header1] : double, eigenvalues
 *  G[header1*header1]  : double, eigenvectors.
 *  	G[i*header+j] is eigenvectors corresponding eigenvalue E[i].
 */
void normal_mode(colloid_base&, const char *, double =-1, bool =true);


/*! density of states core 
 *  @param Gn
 *  	(type)    int
 *  	(input)   number of eigenvalues
 *  @param E
 *  	(type)    double *
 *  	(input)   eigenvalues
 *  @param file
 *  	(type)    const char *
 *  	(input)   filename prefix (can contain subfix) for output
 *
 *  The output is with a subfix ".dos"
 *  This program is called in
 *  	normal_mode
 *  	DOS
 */
void DOS_w(int, double *, const char *);


/*! density of states 
 *  @param ptid
 *  	(type)    colloid_base &
 *  	(input)   [position, ...,  time, id] data
 *  @param file
 *  	(type)    const char *
 *  	(input)   filename prefix (can contain subfix) for output
 *  @param Remove_Drift
 *  	(type)    bool
 *  	(input)   remove drift or not
 *  	(default) remove drift
 *
 *  (input)ptid should contain same number of particles for each frame. 
 *  The program first try to read eigen data from .ev file if found, else
 *  it will call normal_mode to generate everything.
 */
void DOS(colloid_base &, const char *, bool);

/* logarithmic binning of DOS */
void DOSlog(const char *);
void DOS_log(const char *file);

/* DOS obtatined by kernel density estimate */
void DOS_kdensity(const char *file, double h, int nbin);

/* participation ratio and smoothing */
void participation_ratio_w(int, double *, double *, const char *);
void participation_ratio(colloid_base &, const char *, bool);
void participation_ratio_smooth(const char *file, double h);

/* level sapcing distribution and spectral rigidity */
void level_separation_w(int, double *, const char *);
void level_separation(colloid_base &, const char *, bool);
void level_separation_log(const char *);
void level_spacing(const char *file, double h, unfolding_t UNFOLDING);
void level_rigidity(const char *file, double h, unfolding_t UNFOLDING);

/* amplitude distribution of a given mode */
void mode_hist_w(int, double *, int, const char *);
void mode_hist(colloid_base &, int, const char *, bool);

/* obtain covariance matrix */
void covariance(colloid_base &, const char *, bool =true);

/* spatial displacement distribution */
void uur(colloid_base &, const char *, double =-1, bool =true);

/* decomposition of a mode to transverse and longitudinal part */
void transverse_longitudinal(int&, const char *);
void transverse_longitudinal_eps(int&, const char *);

/* plot modes and ouput eps file */
void plot_mode(int&, const char *);
void plot_mode_eps(int&, const char *);

/*! project displacement to normal modes */
void projection(const char *, double =-1, bool =true);


/*! dispersion relation of the system */
void dispersion(const char *);
void plot_dispersion(const char *, int =0);

/*! inverse of dymamic covariance matrix, 
 * the hessian matrix of the system
 */
void K(const char *);
void Kr(const char *);

/*! binary IO of K */
void writeK(int &, double *, const char *);
void readK(int &, double **, const char *);

/* plot K */
void plot_K(const char *file);

/*! potential energy of each frame 
 * should pass the same argument with normal_mode
 */
void U(const char *, double =-1, bool =true);

/* soft spots */
void softspot(const char *file, int Nmode, int Nsmall);

/* spatial correlation of amplitudes of a mode */
void mode_spcorr(int n, const char *file);
void mode_spcorr_all(const char *file, double omegaRange);

/*==================================================================
 * ELLIPSOID NORMAL MODES
 * =================================================================*/
typedef enum { TRANSLATION, ROTATION } freedom_t;
void e_normalmode(const char *file, double a0, double a, double b,
		bool Remove_Drift);

void e_writemp(const int &N, float *mp,const char *file);
void e_readmp(int *Np, float **mp_p, const char *file);

void e_writedcm(int Gn, double *G, const char *file);
void e_readdcm(int *Gnp, double **Gp, const char *file);

void e_writeev(int Gn, double *E, double *G, const char *file);
void e_readev(int *Gnp, double **Ep, double **Gp, const char *file);
void e_readmode(int N, int *Gnp, double *Ep, double **Gp, const char *file);
void e_readE(int *Gnp, double **Ep,const char *file);

void e_plot_mode(const char *file, int n);
void e_plot_mode2(const char *file, double a, double b, double a0, int n);
void e_plot_mode3(const char *file, double a, double b, double a0, int n, double factor);
void e_plot_mode4(const char *file, double a, double b, double a0, int n, 
		freedom_t fdm, double factor);
void e_plot_mode_trans(const char *file, double a, double b, double a0, int n, 
		double factor);
void e_plot_mode_rotat(const char *file, double a, double b, double a0, int n, 
		double factor);
void e_DOS(const char *file);
void e_CDOS(const char *file);
void e_DOS_log(const char *file);
void e_PR(const char *file);
void e_P(const char* file);
void e_corr(const char *file);
void e_plot_config(const char *file, int n, double a, double b,
		bool Tracked);
void e_wrap_angle_normalmode(const char *file, double a0, double a, double b,
		bool Remove_Drift);

void e_softspot(const char *file, double a, double b, double a0, 
		double threshold);
#endif /* NORMAL_MODE_H */
