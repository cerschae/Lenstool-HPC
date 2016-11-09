#include <stdio.h>
 
/* explain function
 * Description:
 *   fixarray::Qsort() is an internal subroutine that implements quick sort.
 *
 * Origin :  
 * 	http://www.d.umn.edu/~gshute/C/examples/quicksort.C.html
 * Return Value: none
 */


void qksort(int ilo, int ihi, double *A) {
    double pivot;		// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    double tempEntry;	// temporary entry used for swapping

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
	if (A[uhi] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo] = A[uhi];
	    A[uhi] = tempEntry;
	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = A[ieq];
		A[ieq] = A[ulo];
		A[ulo] = tempEntry;
		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    qksort(ilo, ieq - 1, A);
    qksort(uhi + 1, ihi, A);
}

void arraySort(double *a,  int len)
{
	qksort(0, len-1, a);
}
