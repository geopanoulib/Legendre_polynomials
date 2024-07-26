/**
 *
 * Copyright (c) 2024, Georgios Panou
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 * Authors: Koci J., Iossifidis C. and Panou G. <geopanou@survey.ntua.gr>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
	register int p, q, t;
	int N, n;
	double theta, pihalf, raddeg, cosnm1x, cosnm2x, cosnx, coeff, *cosn;
	double A, B, C, coss, Pi, pnn, pn0, Tni, *Pn;		
	
	pihalf = 2.0l * atan(1.0l);
	raddeg = pihalf / 90.0l;
	N = 10; /* The maximum degree */
	theta = 30.0; /* The co-latitude in degrees */
	theta *= raddeg;

	cosn = malloc((N + 1) * sizeof(double));
	Pn = malloc((N + 1) * sizeof(double));

	printf("\nN = %d\ntheta = %15.8e deg\n\n%-8s%-16s\n", N, theta / raddeg, "n", "Pn(theta)");
	
	/*
	 * Computes all the multiple angle cosines by using the Chebyshev's method, 
	 * https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Chebyshev_method
     	 *
     	 */
    
	cosnm1x = cos(theta);
	coeff = 2.0l * cosnm1x;
	cosnm2x = 1.0l;
	cosn[0] = 1.0l;
	cosn[1] = cosnm1x;
	for (n = 2; n <= N; n++)
	{
		cosnx = coeff * cosnm1x - cosnm2x;
		cosn[n] = cosnx;
		cosnm2x = cosnm1x;
		cosnm1x = cosnx;
	}
	n = 0;
	Pn[n] = 1.0l;
	
	n = 1;
	Pn[n] = sqrt(3) * cosn[1];

	/* Computes the odd Legendre polynomials (m = 0) */
	pnn = 1.0l;
	coss = cosn[1];
	for (n = 3; n <= N; n += 2)
	{
		p = n - 1;
		A = 1.0l / p;
		B = 1.0l / n;
		coeff = 1.0l - 0.75l * B;
		pnn *= 1.0l - A * coeff;
		Pi = coss;
		q = n + 2;
		t = 3;
		while (p > 0) 
		{
			/* Computes the Tni ratio */
			B = 1.0l - 1.0l / p;
			C = 1.0l + 1.0l / q;
			Tni = B * C;
			Pi *= Tni;
			Pi += cosn[t];
			t += 2;
			p -= 2;
			q += 2;
		}
		Pn[n] = sqrt(2.0l * n + 1) * Pi * pnn;
	}
	
	/* Computes the even Legendre polynomials (m = 0) */
	pn0 = 1.0l;
	pnn = 2.0l;
	coss = cosn[2];
	for (n = 2; n <= N; n += 2)
	{
		p = n - 1;
		/* Computes the Dn_n00 ratio */
		A = 1.0l / p;
		B = 1.0l / n;
		coeff = 1.0l - 0.75l * B;
		pnn *= 1.0l - A * coeff;
		/* Computes the tn ratio */
		A = 1.0l - B;
		pn0 *= A * A;
		Pi = coss;
		q = n + 3;
		t = 4;
		p--;
		while (p > 0)
		{
			/* Computes the Tni ratio */
			B = 1.0l - 1.0l / p;
			C = 1.0l + 1.0l / q;
			Tni = B * C;
			Pi *= Tni;
			Pi += cosn[t];
			t += 2;
			p -= 2;
			q += 2;
		}
		Pn[n] = sqrt(2.0l * n + 1) * (Pi * pnn + pn0);
	}
	/* Printing the Legendre Polynomials */
	for (n = 0; n <= N; n++)
		printf("%-8d%23.16e\n", n, Pn[n]);
	return 0;
}
