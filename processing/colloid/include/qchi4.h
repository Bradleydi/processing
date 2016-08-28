#ifndef QCHI4_H
#define QCHI4_H

#include "colloid_base.h"

/*! four-point correlation suseptibility chi_4(t)
 *  @param pt
 *  	(type) colloid_base&
 *  	untracked position data 
 *  @param a
 *  	(type) double
 *  	overlap function cutoff
 *  @param dim
 *  	(type) int
 *  	dimension
 *
 * Four-point correlation suseptibility is defined in 
 * http://fstarr.web.wesleyan.edu/publications/lssg-jcp.pdf
 * http://www.physics.emory.edu/~weeks/lab/papers/narumi.pdf
 * http://scitation.aip.org/getpdf/servlet/GetPDFServlet?filetype=pdf&id=JCPSA6000112000002000509000001&idtype=cvips&doi=10.1063/1.480541&prog=normal
 * http://spahn.engin.umich.edu/publications/documents/2002_glotzer_58.pdf
 */
void pt_qchi4(colloid_base&, double , int =2);

#endif /* QCHI$_H */
