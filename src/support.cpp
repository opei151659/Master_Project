#include "support.h"

pset_family support_set(pset_family F, pset_family R) {
	pset_family supp = sf_new(NUMOUTPUTS, cube.size);
	int fnum, rnum, var, a, b;
	// Last output word
	// First output word

	unsigned int x, n;
	pset fp, rp, k;

	for (int i = 0; i < NUMOUTPUTS; i++) {
		sf_addset(supp, set_new(cube.size));
		k = supp->data + supp->wsize * (i);
		set_insert(k, NUMINPUTS * 2 + i);
	}



	a = BPI - NUMINPUTS * 2 % BPI;
	b = NUMOUTPUTS - (a + (low - fow - 1) * BPI);

	//printf("%d %d %d %d\n", low, fow, a, b);

	foreachi_set(F, fnum, fp) {
		foreachi_set(R, rnum, rp) {
			//if (output_intersect(fp, rp)) {
			if (OutputIntersection(fp, rp)) {
				var = cdist1(fp, rp);
				if (var == -1) continue;

				x = fp[fow] & rp[fow] & cube.mv_mask[fow];
				if (x) {
					for (int i = BPI - a, n = 0x1 << i; i < BPI; i++, n <<= 1) {
						if (x & n) {
							k = supp->data + supp->wsize * (i - (BPI - a));
							set_insert(k, var * 2);
							set_insert(k, var * 2 + 1);
						}
					}
				}

				if (low != fow) {
					for (int i = fow + 1; i < low; i++) {
						x = fp[i] & rp[i] & cube.mv_mask[i];
						if (x) {
							for (int j = 0, n = 0x1 << j; j < BPI; j++, n <<= 1) {
								if (x & n) {
									k = supp->data + supp->wsize * (a + (i - fow - 1) * BPI + j);
									set_insert(k, var * 2);
									set_insert(k, var * 2 + 1);
								}
							}
						}
					}
					x = fp[low] & rp[low] & cube.mv_mask[low];
					if (x) {
						for (int i = 0, n = 0x1 << i; i < b; i++, n <<= 1) {
							if (x & n) {
								k = supp->data + supp->wsize * (a + (low - fow - 1) * BPI + i);
								set_insert(k, var * 2);
								set_insert(k, var * 2 + 1);
							}
						}
					}
				}
			}
		}
	}

	return supp;
}

/*quick separate against a single output function
using_sup == 1 -> pset & support_set
*/
pcover sep_sup_output(pset_family T, int i, pset sup)
{
	pcover T1;
	pcube p, last, pdest, mask;

	mask = cube.var_mask[cube.output];
	T1 = new_cover(T->count);
	foreach_set(T, last, p) {
		if (is_in_set(p, i)) {
			pdest = GETSET(T1, T1->count++);
			INLINEset_and(pdest, p, sup);
			INLINEset_or(pdest, pdest, mask);
			RESET(pdest, PRIME);
		}
	}
	return T1;
}