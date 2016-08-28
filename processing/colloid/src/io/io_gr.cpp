#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstring>

void write_a0(double &a0, const char *file)
{
	/* test whether a0 exists 
	 * and determine whether should update the information
	 */
	double temp_a0;
	read_a0 (temp_a0, file);

	if ( a0 != temp_a0 )
	{
		char *infofile = getfilename(file, ".info");
		FILE *infof = fopen(infofile, "a");
		FILE_NULL(infof, infofile);

		fprintf(infof, "[ a0 ]\n");
		fprintf(infof, "# lattice constant\n");
		fprintf(infof, "a0 = %f\n\n", a0);

		fclose(infof);

		printf("# lattice constant a0 = %f has been written to file \"%s\".\n",
				a0, infofile);

		free(infofile);
	}
}


void read_a0(double &a0, const char *file)
{
	unsigned char is_default=1;
	/*================================
	 * check .info file exists or not?
	 *================================*/

	char *infofile = getfilename(file, ".info");
	FILE *infof = fopen(infofile, "r");
	if ( infof != NULL )
	{
		const int buffersize=80;
		char *buffer=(char *)malloc(buffersize*sizeof(char));
		POINTER_NULL(buffer);

		const char *section="[ a0 ]\n";
		char *key=(char *)malloc(buffersize*sizeof(char));
		POINTER_NULL(key);
		float tmp_a0;
		
		unsigned char get_info=0;

		while ( get_info ==0 && fgets(buffer, buffersize, infof) != NULL )
		{
			//printf("# %s\n", buffer);
			// compare with string section="[ a0 ]\n"
			if ( strcmp(buffer, section) == 0 ) // equal
			{
				get_info=1;
				// read parameter
				while ( fgets(buffer, buffersize, infof) != NULL )
				{
					if ( buffer[0] != '#' ) // not a comment line
					{
						if ( buffer[0] != '\n' ) // not end of section
						{
							// printf ("# getting a0...\n");
							if ( sscanf(buffer, "%s = %f", key, &tmp_a0) == 2 )
							{
								a0=(double)tmp_a0;
								printf("# File \"%s\" loaded.\n"
									   "#     a0 = %f\n", infofile, a0);
								is_default=0;
								break;
							}
						}
						else if ( buffer[0] == '\n' ) // end of section
							break;
					}
				}
			}
		}
		fclose(infof);
		
		free(buffer);
		free(key);
	}
	
	if ( is_default )
	{
		a0=-1;
		printf("# File \"%s\" doesn't exist.\n# set a0 = -1\n"
				, infofile);
	}

	free(infofile);
}
