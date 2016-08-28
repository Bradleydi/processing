#include <miscellaneous.h>

/* copy from numerical recipes in C++ version 2.10 */
/* WARNING 1: array is also changed 
 * WARNING 2: No check for k compared with 0 and n*/
#define SWAP(a, b) temp=(a);(a)=(b);(b)=temp;

/* Returns the kth smallest value in the array arr[1..n]. The input array will 
 * be rearranged to have this value in location arr[k], with all smaller 
 * elements moved to arr[1..k-1] (in arbitrary order) and all larger elements 
 * in arr[k+1..n] (also in arbitrary order)
 */

double select_NR(int n, double *arr, int k)
{
	int i, ir, j, l, mid;
	double a, temp;

	l=0;
	ir=n-1;
	--k; // kth is actuall in C the (k-1)th

	for (;;) {
		if (ir <= l+1) {		//Active partition contains 1 or 2 elements.
			if (ir == l+1 && arr[ir] < arr[l]) {		// Case of 2 elements.
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1; // Choose median of left, center, and right el-
			SWAP(arr[mid],arr[l+1]) // ements as partitioning element a. Also
			if (arr[l] > arr[ir]) { //rearrange so that arr[l] ≤ arr[l+1],
				SWAP(arr[l],arr[ir])  // arr[ir] ≥ arr[l+1].
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
			}
			i=l+1; // Initialize pointers for partitioning.
			j=ir;
			a=arr[l+1]; // Partitioning element.
			for (;;) {  // Beginning of innermost loop
				do i++; while (arr[i] < a); // Scan up to find element > a.
				do j--; while (arr[j] > a); // Scan down to find element < a.
				if (j < i) break;  // Pointers crossed. Partitioning complete.
				SWAP(arr[i],arr[j])
				/*
				do i++; while (arr[i] < a); // Scan up to find element > a.
				do j--; while (arr[j] > a); // Scan down to find element < a.
				if (j < i) break;  // Pointers crossed. Partitioning complete.
				SWAP(arr[i],arr[j])
				// Here do j--; while (arr[j] > a); has a bug, it will
				// make j<0, so that arr[j] is wrong
				*/
			}
			arr[l+1]=arr[j]; // Insert partitioning element.
			arr[j]=a;
			if (j >= k) ir=j-1; // Keep active the partition that contains the
			if (j <= k) l=i; // kth element.
		}
	}
}

#undef SWAP

double select(int n, double *arr, int k)
{
	double *copy = Malloc(double, n); POINTER_NULL(copy);
	int i;
	for (i=0; i<n; i++) copy[i] = arr[i];
	double kth = select_NR(n, copy, k);
	free(copy);
	return kth;
}

/* select with given copy array allocated, this is useful when we need
 * to use select many times with array same length.
 * Hence in this case, we don't need to allocated the copy array again and again
 */
double select_wcp(int n, double *arr, double *copy, int k)
{
	int i;
	for (i=0; i<n; i++) copy[i] = arr[i];
	return select_NR(n, copy, k);
}
