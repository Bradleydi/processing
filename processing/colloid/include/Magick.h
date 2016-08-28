#ifndef MAGICK_H
#define MAGICK_H

#include "colloid_base.h"


/*! A function to read image file
 *
 *	@param cb
 *		store the read data, should be initialized by colloid_base cb;
 *	@param filename
 *		name of the image file to input
 *  @param not_in_use
 *      true: cb should be reserved memory to use
 *      flase: cb has been reserved memory
 *
 *  Caution: not_in_use parameter is used to avoid reallocate memory for cb,
 *           cb should be only allocate memory once.
 *           In the program, it's also assumed that the allocated memory is
 *           same size as the image to input. No error message will be showed
 *           when dismatched. So the users should be very careful when dealing
 *           with this function.
 *
 */
void readimage(colloid_base &, 
		       const char  *,
			   const bool & =true);

/*! Show xy data in an image
 *  @param p
 *  	a data containing xy position information, assuming as the first and
 *  	second column
 *	@param factor
 *		a factor to resize the image to make the image better
 *	@param filename
 *		name of the image file to output
 *
 */
/*
void showposition(colloid_base &, 
	              const int =1,
		          const char * ="tmp.tiff");
*/
#endif /* MAGICK_H */
