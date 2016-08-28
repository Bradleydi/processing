#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

char * getfilename(const char *file, const char *subfix)
{
	int count=0, sc=0, i;
	while (file[count]!='\0')
		++count;
	i=count-1;
	while (file[i]!='.' && i>0)
		--i;
	if (i!=0)
		count=i;

	while (subfix[sc]!='\0')
		++sc;
	
	char *filename=(char *)malloc((count+sc+1)*sizeof(char));
	POINTER_NULL(filename);
	for (i=0; i<count; i++)
		filename[i]=file[i];
	i=0;
	while (subfix[i]!='\0')
	{
		filename[count+i]=subfix[i];
		++i;
	}
	filename[count+i]='\0';
	return filename;
}
