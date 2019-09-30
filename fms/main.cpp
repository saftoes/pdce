#include "arb.h"
#include "arb_poly.h"

#include <sys/time.h>
#include <gmp.h>
#include <unistd.h>
#include <getopt.h>

#define PREC 	300
#define DIGIT	20

#define _DYNAMIC_BOUND
#define _SLIENT

// Computing v(j)
mp_limb_t v(mp_limb_t j) {

	mp_bitcnt_t nbits = mpn_popcount(&j, 1);
	if (nbits%2 == 0) return 1;
	else return -1;
}

// Computing p_{n,z}
// If macro _DYNAMIC_BOUND used, use dynamic truncation error bound paratermized by delta
void p(arb_t res,
	mp_limb_t n, mp_limb_t z, mp_limb_t delta) {

	ulong loop = n_pow(2, z);
	arb_zero(res);

	arb_t temp, watch, threshold, den;

	arb_init(temp);
	arb_init(watch);
	arb_init(threshold);
	arb_init(den);

	arb_ui_pow_ui(threshold, 2, delta, PREC);
	arb_ui_div(threshold, 1, threshold, PREC);
	arb_set_ui(den, 2*loop);

	ulong j = 0;
	for (; j < loop; j+=2) {
		arb_zero(watch);

		arb_ui_div(temp, 2*j + 1, den, PREC);
		arb_sub_ui(temp, temp, 1, PREC);
		arb_neg(temp, temp);					// Not necessary if n is known to be even
		arb_pow_ui(temp, temp, n, PREC);
		arb_addmul_si(watch, temp, v(2*j + 1), PREC);

		if (loop != 1) {
			arb_ui_div(temp, 2*j + 3, den, PREC);
			arb_sub_ui(temp, temp, 1, PREC);
			arb_neg(temp, temp);				// Not necessary if n is known to be even
			arb_pow_ui(temp, temp, n, PREC);
			arb_addmul_si(watch, temp, v(2*j + 3), PREC);
		}		

#ifdef _DYNAMIC_BOUND
		arb_mul_ui(temp, watch, (loop - j + 1)/2, PREC);
		arb_abs(temp, temp);
		if (arb_le(temp, threshold) == 1) {		// Truncation error bound is less than threshold
#ifndef _SLIENT
			arb_printd(watch, DIGIT);		printf("\n");
			arb_printd(temp, DIGIT);		printf("\n");
			arb_printd(threshold, DIGIT);	printf("\n");
#endif
			break;
		}
		else arb_sub(res, res, watch, PREC);
#else
		arb_sub(res, res, watch, PREC);
#endif
	}
#ifndef _SLIENT
	printf("End with j = %ld (< %ld?)\n", j, loop);
#endif

	arb_clear(temp);
	arb_clear(watch);
	arb_clear(threshold);
	arb_clear(den);
}
/*
void fast_arb_poly_pow_ui(arb_poly_t res,
	arb_poly_t poly, mp_limb_t exp, slong prec) {

	arb_poly_one(res);

	int counter = 1;
	struct timeval start, end;

	for (mp_limb_t e = exp; e > 0; e >>= 1) {

		printf("%d\n", counter++);
		gettimeofday(&start, NULL);

		if ((e & 0x1) == 0) {
			arb_poly_pow_ui(poly, poly, 2, prec);

			gettimeofday(&end, NULL);
			printf("    Time used: %f ms\n", 
				(end.tv_sec-start.tv_sec)*1000.0 + (end.tv_usec-start.tv_usec)/1000.0);

			continue;
		}

		arb_poly_mul(res, res, poly, prec);

		gettimeofday(&end, NULL);
		printf("    Time used: %f ms\n", 
			(end.tv_sec-start.tv_sec)*1000.0 + (end.tv_usec-start.tv_usec)/1000.0);
		gettimeofday(&start, NULL);

		arb_poly_pow_ui(poly, poly, 2, prec);

		gettimeofday(&end, NULL);
		printf("    Time used: %f ms\n", 
			(end.tv_sec-start.tv_sec)*1000.0 + (end.tv_usec-start.tv_usec)/1000.0);

	}
}
*/

// Computing probability generating function of *big Z*
void pgf(arb_poly_t poly,
	mp_limb_t n, mp_limb_t w, mp_limb_t m, mp_limb_t delta) {

	arb_t coeff;
	arb_init(coeff);

	arb_poly_t pgfz;
	arb_poly_init(pgfz);
	arb_poly_fit_length(pgfz, w);

	for (ulong z = 0; z < w; ++z) {
		p(coeff, n, z, delta);
		arb_poly_set_coeff_arb(pgfz, z, coeff);
	}

	struct timeval start, end;
	gettimeofday(&start, NULL);

	arb_poly_pow_ui(poly, pgfz, m, PREC);
	//fast_arb_poly_pow_ui(poly, pgfz, m, PREC);

	gettimeofday(&end, NULL);
	printf("Powering time: %f ms\n", 
		(end.tv_sec-start.tv_sec)*1000.0 + (end.tv_usec-start.tv_usec)/1000.0);

	arb_clear(coeff);
	arb_poly_clear(pgfz);
}

// Computing probability mass function of *big Z*
void pmf(arb_t res, arb_poly_t poly, mp_limb_t Z) {

	arb_poly_get_coeff_arb(res, poly, Z);
}

// Find the smallest N0 satisfying DP
mp_limb_t find_N0(mp_limb_t w, mp_limb_t m, double epsilon, mp_limb_t delta, mp_limb_t init) {

	arb_t _epsilon, _epsilon_inv, _delta, temp1, temp2, left, right;
	arb_poly_t D1, D2;

	arb_init(_epsilon);
	arb_init(_epsilon_inv);
	arb_init(_delta);
	arb_init(temp1);
	arb_init(temp2);
	arb_init(left);
	arb_init(right);
	arb_poly_init(D1);
	arb_poly_init(D2);

	arb_set_d(temp2, epsilon);
	arb_exp_invexp(_epsilon, _epsilon_inv, temp2, PREC);	// _epsilon = e^{epsilon}
															// _epsilon_inv = e^{-epsilon}
	arb_set_si(temp2, -delta);
	arb_set_ui(_delta, 2);
	arb_pow(_delta, _delta, temp2, PREC);					// _delta = 2^{-delta}

	ulong floor = 1, ceil = init;
	mp_limb_t count = 0;
	while(count != ((w - 1)*m + 1)) {
		// Probability generation function of D1 where |D1| = ceil
		pgf(D1, ceil, w, m, PREC);		printf("D1 done\n");
		// Probability generation function of D2 where |D2| = ceil + 1
		pgf(D2, ceil + 1, w, m, PREC);	printf("D2 done\n");

		count = 0;
		for (ulong k = 0; k <= (ulong)((w - 1)*m); ++k) {
			pmf(temp1, D1, k);	pmf(temp2, D2, k);		// temp1 = Pr[Z_N = K], temp2 = Pr[Z_{N+1} = K]

			arb_mul(left, _epsilon_inv, temp2, PREC);
			arb_sub(left, left, _delta, PREC);			// left = _epsilon_inv*Pr[Z_{N+1} = K] - _delta

			arb_mul(right, _epsilon, temp2, PREC);
			arb_add(right, right, _delta, PREC);		// right = _epsilon*Pr[Z_{N+1} = K] + _delta

			if (arb_ge(temp1, left) && arb_le(temp1, right)) count++;
			else {
				floor = ceil;
				ceil *= 2;
				break;
			}
		}
		printf("Count: %ld\n", count);
	}

	printf("Suspected interval: [%ld, %ld]\n\n", floor, ceil);

	bool flag;
	ulong pivot;

	for ( ; floor != ceil - 1; ) {
		pivot = (floor + ceil)/2;
		printf("pivot = %ld\n", pivot);

		// Probability generation function of D1 where |D1| = pivot
		pgf(D1, pivot, w, m, PREC);			printf("D1 done\n");
		// Probability generation function of D2 where |D2| = pivot + 1
		pgf(D2, pivot + 1, w, m, PREC); 	printf("D2 done\n");

		flag = true;
		for (ulong k = 0; k <= (ulong)((w - 1)*m); ++k) {
			pmf(temp1, D1, k);	pmf(temp2, D2, k);		// temp1 = Pr[Z_N = K], temp2 = Pr[Z_{N+1} = K]

			arb_mul(left, _epsilon_inv, temp2, PREC);
			arb_sub(left, left, _delta, PREC);			// left = _epsilon_inv*Pr[Z_{N+1} = K] - _delta

			arb_mul(right, _epsilon, temp2, PREC);
			arb_add(right, right, _delta, PREC);		// right = _epsilon*Pr[Z_{N+1} = K] + _delta

			if (!(arb_ge(temp1, left) && arb_le(temp1, right))) {	// Pivot does not satisfy DP
				flag = false;
				floor = pivot;
				break;
			}
		}

		if (flag) ceil = pivot;

		printf("floor = %ld, ceil = %ld\n\n", floor, ceil);
	}

	arb_clear(_epsilon);
	arb_clear(_epsilon_inv);
	arb_clear(_delta);
	arb_clear(temp1);
	arb_clear(temp2);
	arb_clear(left);
	arb_clear(right);
	arb_poly_clear(D1);
	arb_poly_clear(D2);

	return ceil;
}

/*
// Debug subroutines
void _check() {

	arb_t check1, check2, temp;
	arb_init(check1);
	arb_init(check2);
	arb_init(temp);

	// Check p_{n,z}
	arb_zero(check1);
	arb_zero(check2);
	for (ulong j = 0; j < 32; ++j) {
		p(temp, 10000, j, PREC);
		arb_add(check1, check1, temp, PREC);
		arb_addmul_ui(check2, temp, j, PREC);
	}

	arb_set_d(temp, 0.77351*10000);
	arb_log_base_ui(temp, temp, 2, PREC);

	printf("\n[Sum] 1 ~= ");	arb_printd(check1, DIGIT);
	printf("\n[Exp] ");			arb_printd(temp, DIGIT);
	printf(" ~= ");				arb_printd(check2, DIGIT);
	printf("\n");

	arb_clear(check1);
	arb_clear(check2);
	arb_clear(temp);
}
*/

int main(int argc, char *argv[])
{
	// Parse command line args
	struct option opts[] = {
		// Three modes
		{"individual", 	0, NULL, 'i'},
		{"aggregation", 0, NULL, 'a'},
		{"search", 		0, NULL, 's'},
		// Program parameters
		{"n", 		1, NULL, 'n'},
		{"w", 		1, NULL, 'w'},
		{"m", 		1, NULL, 'm'},
		{"epsilon", 1, NULL, 'e'},
		{"delta",	1, NULL, 'd'},
		// Other parameters
		{"out",	 	1, NULL, 'o'},
		{"ceil", 	1, NULL, 'c'},
		{"help", 	0, NULL, 'h'}
	};

	int op;
	bool individual = false, aggregation = false, search = false;
	int n = -1, w = -1, m = -1;
	double epsilon = -1;
	int delta = 40;
	int ceil = 8192;
	FILE *fp = NULL;

	while((op = getopt_long(argc, argv, "iasn:w:m:e:d:o:c:h", opts, NULL)) != -1){
		switch(op){
			case 'i': 	individual = true; 				break;
			case 'a': 	aggregation = true; 			break;
			case 's':	search = true;					break;

			case 'n':	n = atoi(optarg);				break;
			case 'w':	w = atoi(optarg);				break;
			case 'm':	m = atoi(optarg);				break;
			case 'e':	epsilon = atof(optarg);			break;
			case 'd':	delta = atoi(optarg);			break;

			case 'o':	fp = fopen(optarg, "w");		break;
			case 'c':	ceil = atoi(optarg);			break;
			case '?':	printf("Invalid option!\n");	break;
			case 'h':
			default:	return -1;
		}
	}

	if (!individual && !aggregation && !search) {
		printf("Incomplete command line args! Please specify at least one task: -i/-a/-s\n");
		return -1;
	}

	if (fp == NULL) fp = fopen("fms.log", "w");
	arb_t res;
	arb_poly_t poly;

	if (individual) {
		if (n <= 0 || w <= 0) {
			printf("Invalid parameters for -i. Proper -n and -w are required.\n");
			return -1;
		}

		arb_init(res);
		arb_poly_init(poly);

		fprintf(fp, "Computing probability generating function of z\n");
		fprintf(fp, "n = %d, w = %d\n", n, w);
		fprintf(fp, "--------------------------------------------------------\n");

		for (ulong z = 0; z < (ulong)w; ++z) {
			p(res, n, z, PREC);
			fprintf(fp, "z = %ld: ", z);
			arb_fprintd(fp, res, DIGIT);
			fprintf(fp, "\n");
		}

		fclose(fp);
		arb_clear(res);
		arb_poly_clear(poly);
		return 0;
	}

	if (aggregation) {
		if (n <= 0 || w <= 0 || m <= 0) {
			printf("Invalid parameters for -a. Proper -n, -w and -m are required.\n");
			return -1;
		}

		arb_init(res);
		arb_poly_init(poly);

		pgf(poly, n, w, m, PREC);
		fprintf(fp, "Computing probability mass function of Z\n");
		fprintf(fp, "n = %d, w = %d, m = %d\n", n, w, m);
		fprintf(fp, "Degree of PGF of Z: %ld\n", arb_poly_degree(poly));
		fprintf(fp, "--------------------------------------------------------\n");
		
		//for (ulong i = 0; i <= (ulong)arb_poly_degree(poly); ++i) {
		for (ulong i = 0; i <= (ulong)((w - 1)*m); ++i) {
			pmf(res, poly, i);
			fprintf(fp, "Z = %ld: ", i);
			arb_fprintd(fp, res, DIGIT);
			fprintf(fp, "\n");
		}

		fclose(fp);
		arb_clear(res);
		arb_poly_clear(poly);
		return 0;
	}

	if (search) {
		if (w <= 0 || m <= 0 || epsilon <= 0 || delta <= 0 || ceil <= 0) {
			printf("Invalid parameters for -s. Proper -w, -m, -e are required (-d and -c are optional).\n");
			return -1;
		}

		ulong N0 = find_N0(w, m, epsilon, delta, ceil);

		fprintf(fp, "w = %d, m = %d, epsilon = %f, delta = 2^{-%d}\n", w, m, epsilon, delta);
		fprintf(fp, "Set initial ceiling to %d\n", ceil);
		fprintf(fp, "N0 = %ld\n", N0);

		printf("w = %d, m = %d, epsilon = %f, delta = 2^{-%d}\n", w, m, epsilon, delta);
		printf("Set initial ceiling to %d\n", ceil);
		printf("N0 = %ld\n", N0);

		fclose(fp);
		return 0;
	}
}