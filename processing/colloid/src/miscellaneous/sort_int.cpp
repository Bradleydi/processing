#include "miscellaneous.h"

using namespace std;

void sort_int(int n, int *p)
{
	int *temp=(int *)malloc(n*sizeof(int));
	POINTER_NULL(temp);

	Merge_sort_int(p, temp, 0, n-1);
	free(temp);
}


// sort array p from index l to index r
void Merge_sort_int(int *p, int *temp, int l, int r)
{
	int i, j, k, m;
	
	if ( r-l > 0 )
	{
		// divide
		m=(r+l)/2;
		Merge_sort_int(p, temp, l, m);
		Merge_sort_int(p, temp, m+1, r);

		//copy
		for (i=l; i<=m; i++)
			temp[i]=p[i];
		for (j=m+1; j<=r; j++)
			temp[j]=p[j];

		// merge
		i=l;
		j=m+1;
		for (k=l; k<=r; k++)
		{
			if ( temp[i] < temp[j] )
			{
				p[k]=temp[i];
				if (i==m) // last i, so all remaining are j's
				{
					for (i=k+1; i<=r; i++)
						p[i]=temp[j++];
					break;
				}
				else
					++i;
			}
			else
			{
				p[k]=temp[j];
				if (j==r) // last i, so all remaining are j's
				{
					for (j=k+1; j<=r; j++)
						p[j]=temp[i++];
					break;
				}
				else
					++j;
			}
		}
	}
}
