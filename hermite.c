#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "makespl.h"
#include "gaus/piv_ge_solver.h"

double hermite(int n, double x) {
	if (n == 0) {
		return 1;
	}
	if (n == 1) {
		return 2*x;
	}
	return 2 * x * hermite(n - 1,x) - 2 * (n - 1) * hermite( n - 2,x); //przesuwamy indeksy h(x) o jeden w lewo, stąd n - 1 ; n - 2
}

double d1(int n, double x) {
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		return 2;
	}
	return 2 * hermite(n -1,x) + 2 * x * d1(n -1,x) - 2 * (n -1) * d1(n - 2,x);
}

double d2(int n, double x) {
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		return 0;
	}
	return 4 * d1(n -1,x) + 2 * x * d2(n -1,x) - 2 * (n-1) * d2(n - 2,x);
}

double d3(int n, double x) {
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		return 0;
	}
	return 6 * d2(n -1,x) + 2 * x * d3(n -1,x) - 2 * (n -1) * d3(n - 2,x);
}


void make_spl(points_t* pts, spline_t* spl) {
	matrix_t* eqs = NULL;
	double* x = pts->x;
	double* y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1]; 
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
	char* nbEnv = getenv("APPROX_BASE_SIZE");

	if (nbEnv != NULL && atoi(nbEnv) > 0)
		nb = atoi(nbEnv);

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, hermite(i, x[k]) * hermite(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * hermite(j, x[k]));
	}

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i * (b - a) / (spl->n - 1);
			xx += 10.0 * DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i] += ck * hermite(k, xx);
				spl->f1[i] += ck * d1(k, xx);
				spl->f2[i] += ck * d2(k, xx);
				spl->f3[i] += ck * d3(k, xx);
			}
		}
	}
}
