#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

void pformat();

int* getframeindex(const char *string, int lower, int upper)
{

	if (lower>upper)
		pERROR("wrong upper & lower bound");

	int *frameindex=Malloc(int, 3); POINTER_NULL(frameindex);

	//==============================================================
	// check the string, whether is composed by integers and '+' '-' ':' only
	int i, N=0;
	while (string[N]!='\0')
		++N;

	int Ncolon=0;
	for (i=0; i<N; i++)
	{
		if(string[i]==':')
			++Ncolon;
		else if (string[i]<'0')
		{
			if (string[i]!='+' && string[i]!='-')
				pformat();
		}
		else if (string[i]>'9')
			pformat();
	}
	//==============================================================
	char *str=Malloc(char, N+1); POINTER_NULL(str);
	for (i=0; i<=N; i++)
		str[i]=string[i];
	
	if (Ncolon==0) // no colon, so it's an integer only
	{
		frameindex[0]=atoi(string);
		frameindex[1]=frameindex[0];
		frameindex[2]=1;
	}
	else if (Ncolon==1) // 1 colon, so 2 integer, lower and upper limit
	{
		int pcolon=0;
		while (string[pcolon]!=':')
			++pcolon;

		str[pcolon]='\0';
		
		if ( pcolon==0 ) // no lower limit specified, meaning from smallest
			frameindex[0]=lower;
		else
			frameindex[0]=atoi(str);

		if ( pcolon==N-1 ) // no upper limit specified, meaning to greatest
			frameindex[1]=upper;
		else
			frameindex[1]=atoi(str+pcolon+1);

		frameindex[2]=1;
	}
	else if (Ncolon==2) // 2colons, so 3 integers, lower, increment, upper 
	{
		int pcolon1=0;
		while (string[pcolon1]!=':')
			++pcolon1;
		int pcolon2=pcolon1+1;
		while (string[pcolon2]!=':')
			++pcolon2;

		str[pcolon1]='\0';
		str[pcolon2]='\0';
		
		if ( pcolon1==0 ) // no lower limit specified, meaning from smallest
			frameindex[0]=lower;
		else
			frameindex[0]=atoi(str);

		if (pcolon1+1==pcolon2) // no increment specified, default to 1
			frameindex[2]=1;
		else
			frameindex[2]=atoi(str+pcolon1+1);

		if (pcolon2==N-1) // no upper lime specified
			frameindex[1]=upper;
		else
			frameindex[1]=atoi(str+pcolon2+1);
	}
	else
		pformat();

	free(str);

	if (frameindex[0]<lower)
	{
		printf("# Reset the lower bound to %d\n", lower);
		frameindex[0]=lower;
	}
	if (frameindex[1]>upper)
	{
		printf("# Reset the upper bound to %d\n", upper);
		frameindex[1]=upper;
	}
	if (frameindex[2]>(upper-lower))
	{
		printf("# Reset the increment to %d\n", upper-lower);
		frameindex[2]=upper-lower;
	}

	// check
	if (frameindex[0]>frameindex[1])
	{
		printf("# Error: Incorrect index, please reset it:\n"
			   "#        Lower bound = %d\n"
			   "#        Upper bound = %d\n", lower, upper);
		exit (1);
	}
	if (frameindex[2]<=0)
		pERROR("Incorrect increment.");

	return frameindex;
}


void pformat()
{
	printf("string 'frameindex' format:\n"
		   "\tN\n"
		   "\t\tset chosen frame as N (only one frame)\n"
		   "\tN1:N2\n"
		   "\t\tN1: lower bound, if not specified, set to lower\n"
		   "\t\tN2: upper bound, if not specified, set to upper\n"
		   "\t\tThe increment is set to 1\n"
		   "\t\tA ':' means choosing all frames.\n"
		   "\tN1:c:N2\n"
		   "\t\tN1: lower bound, if not specified, set to lower\n"
		   "\t\tN2: upper bound, if not specified, set to upper\n"
		   "\t\tc: increment, if not specified, set to 1\n");
	exit(1);
}
