#include "normalmode.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void K(const char *file)
{
	int Gn;
	double *G;
	readdcm(Gn, &G, file);

	dpotri(G, Gn);
	writeK(Gn, G, file);
}
