/*
  Author:  Pate Williams (c) 1997 & James Wanless (c) 2009-10

  Algorithm 10.3.3 (Lenstra's ECM). See "A Course
  in Computational Algebraic Number Theory" by
  Henri Cohen page 488.
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <gmpxx.h>
#include <gmp.h>

using namespace std;

struct point {mpz_class x, y;};

int Rabin_Miller(mpz_class n)
/* given an integer n >= 3 returns 0 composite, 1 probably prime */
{
  return mpz_probab_prime_p(n.get_mpz_t(), 10);
}

mpz_class gcd(mpz_class a, mpz_class b)
/* returns greatest common divisor of a and b */
{
  mpz_t temp;
  mpz_init(temp);
  mpz_gcd(temp, a.get_mpz_t(), b.get_mpz_t());
  mpz_class temp_class(temp);
  mpz_clear(temp);
  return temp_class;
}

mpz_class nextp(mpz_class n)
/* returns next prime after n */
{
	mpz_t temp;
	mpz_init(temp);
	mpz_nextprime(temp, n.get_mpz_t());
	mpz_class temp_class(temp);
	mpz_clear(temp);
	return temp_class;
}

mpz_class inverse(mpz_class a, mpz_class b)
/* returns inverse of a modulo b or 0 if it does not exist */
{
	mpz_t temp;
	mpz_init(temp);
	if (!mpz_invert(temp, a.get_mpz_t(), b.get_mpz_t()))
		mpz_set_si(temp, 0);
	mpz_class temp_class(temp);
	mpz_clear(temp);
	return temp_class;
}

mpz_class modpos(mpz_class a, mpz_class b)
/* returns a modulo b, strictly non-negative */
{
	mpz_class temp_class;
	temp_class = a % b;
	if (temp_class < 0)
		temp_class += b;
	return temp_class;
}

int addition_1(mpz_class n, struct point P1, struct point P2, struct point *P3)
/* affine coords */
/* returns 1 if P1 = -P2 therefore P1 + P2 = O, 0 otherwise */
/* P1 != P2 */
{
	mpz_class delta_x;
	mpz_class delta_y;
	mpz_class m;
	
	delta_x = modpos(P2.x - P1.x, n);
	delta_y = modpos(P2.y - P1.y, n);

	if (P1.x == P2.x && ((P1.y + P2.y) == 0 || (P1.y + P2.y) == n || (P1.y + P2.y) == 2 * n)) {
		P3->x = 0, P3->y = 1;
	return 1;
	}
	
	/* calculate m = (y2 - y1)(x2 - x1) ^ -1 mod n */
	m = modpos(delta_y * inverse(delta_x, n), n);

	/* calculate x3 = m ^ 2 - (x1 + x2) mod n */
	P3->x = modpos(m * m - (P1.x + P2.x), n);

	/* calculate y3 = m(x1 - x3) - y1 mod n */
	P3->y = modpos(m * (P1.x - P3->x) - P1.y, n);
	
	return 0;
}

void addition_2(mpz_class a, mpz_class n, struct point P1, struct point *P3)
/* affine coords */
/* P1 == P2 */
{
	mpz_class m;

	/* calculate m = (3x1 ^ 2 + a)(2y1) ^ -1 mod n */
	m = modpos((3 * P1.x * P1.x + a) * inverse(2 * P1.y, n), n);

	/* calculate x3 = m ^ 2 - 2x1 mod n */
	P3->x = modpos(m * m - 2 * P1.x, n);

	/* calculate y3 = m(x1 - x3) - y1 mod n */
	P3->y = modpos(m * (P1.x - P3->x) - P1.y, n);
}

int multiply(mpz_class a, mpz_class k, mpz_class n, struct point P, struct point *R, mpz_class *d)
/* binary ladder */
/* returns -1 if O encountered, 0 if divisor not found, 1 otherwise */
{
	int value = 1;
	struct point A, B, C;

	/*  A = P */
	A = P;
	/* B = O = (0, 1) the point at infinity */
	B.x = 0, B.y = 1;

	while (value && k > 0) {
		if (k % 2 != 0) {
			*d = gcd(modpos(B.x - A.x, n), n);

			k--;
			value = (*d == 1 || *d == n);
			
			if (A.x == 0 && A.y == 1);
			else if (B.x == 0 && B.y == 1) B = A;
			else if (value) {
				addition_1(n, A, B, &C);
				B = C;
			}
		}
		else {
			*d = gcd(modpos(2 * A.y, n), n);

			k >>= 1;
			value = (*d == 1 || *d == n);
			
			if (value) {
				addition_2(a, n, A, &C);
				A = C;
			}
		}
	}
	
	*R = B;
	R->x = modpos(R->x, n);
	R->y = modpos(R->y, n);

	if (R->x == 0 && R->y == 1) return -1;
	
	return !value;
}


int LenstrasECM(mpz_class *zN, mpz_class *zg)
{
	int found = 0;
	mpz_class B = 10000l;
	mpz_class l, q, q1, newq;
	struct point x, y;
	mpz_class za, zd = 0;
	long curve = 0;

	do {
		
		do {
			za = modpos(rand(), *zN);
			x.x = 0;
			x.y = 1;
	
			curve++;
			cout << "B=" << B << ", curve#" << curve << ", za=" << za << "          \r";
			fflush(stdout);

			newq = 2;
			for (; newq < B && found != 1;) {
				q = newq;
				q1 = q;
				l = B / q;
				while (q1 <= l)
					q1 *= q;

				found = multiply(za, q1, *zN, x, &y, &zd);

				x.x = y.x;
				x.y = y.y;
//  cout << "X=" << x.zx << "," << x.zy << "\n";
				
				newq = nextp(q);
			}

			*zg = gcd(zd, *zN);
		
		} while (curve < 50 && (*zg == *zN || *zg == 1));
  
		B = B * 10;
		curve = 0;

	} while (*zg == *zN || *zg == 1);
  
	cout << "\n";
  
	return 0;
}

int main(void)
{
	char answer[256];
	mpz_class zN = 0, zg = 0;

	do {
		printf("enter the number to be factored below:\n");
		cin >> zN;
	  
		if (Rabin_Miller(zN))
			printf("number is prime\n");
		else {
			printf("factors:\n");
			while (zN % 2 == 0) {
				cout << "2\n";
				zN /= 2;
			}
			while (zN % 3 == 0) {
				cout << "3\n";
				zN /= 3;
			}
			while (!Rabin_Miller(zN)) {
				LenstrasECM(&zN, &zg);
				cout << zg << "\n";
				zN /= zg;
			}
			cout << zN << "\n";
		}
		
		printf("\nanother number (n or y)? \n");
		scanf("%s", answer);
	
	} while (tolower(answer[0]) == 'y');
	
	return 0;
}
