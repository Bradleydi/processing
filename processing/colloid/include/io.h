#ifndef IO_H
#define IO_H

#include "colloid_base.h"

/*! show file size
 * @param filesize
 * 		file size
 */
void show_file_size(unsigned int);


/*! A function to read inputted gdf data
 *
 * GDF data is used in colloid community, which is actually a binary data.
 *	@param cb
 *		store the read data, initialized by colloid_base cb;
 *	@param filename
 *		name of the file to input
 *
 */
void readgdf(colloid_base &, const char *);


/*! A function to write data to .gdf file
 *
 * GDF data is used in colloid community, which is actually a binary data.
 *	@param cb
 *		store the data to write
 *	@param filename
 *		name of the file to output
 *
 */
void writegdf(colloid_base &, const char *);


void write_a0 (double &, const char *);
void read_a0 (double &, const char *);

#endif /* IO_H */
