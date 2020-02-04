#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* funció que resol un sistema d'equacions amb una matriu triangular inferior passada per paràmetre i un vector b també donat*/
int resLinf(int n, double **L, double *b, double *x, double tol){
	int i, k;
	double sum, numerator, denominator;
	numerator = b[0];
	if (fabs(numerator) < tol) numerator = 0.;
	denominator = L[0][0];
	/* si el denominador és nul, el sistema és incompatible i no pot ser resolt*/
	if (fabs(denominator) < tol) return 1;
	x[0] = numerator / denominator;
	for (i = 1; i < n; i++){
		sum = 0.;
		for (k = 0; k < i; k++){
			sum += L[i][k] * x[k];
		}
		numerator = (b[i] - sum);
		if (fabs(numerator) < tol) numerator = 0.;
		denominator = L[i][i];
		/* si el denominador és nul, el sistema és incompatible i no pot ser resolt*/
		if (fabs(denominator) < tol) return 1;
		x[i] = numerator / denominator;
		if (fabs(x[i]) < tol) x[i] = 0.;
	}
	return 0;
}

/*funció que rep una matriu i un vector dels quals obté el producte, que guarda com un altre vector*/
void prodMatVect(int m, int n, double ** A, double * x, double * y){
	int i, j;
	double sum;
	for (i = 0; i < m; i++){
		sum = 0.;
		for (j = 0; j < n; j++){
			sum += A[i][j] * x[j];
		}
		y[i] = sum;
	}

	printf("El vector producte resulta ser:\n");
	for (i = 0; i < m; i++){
		printf("%e\n", y[i]);
	}
}


/* funció que resol un sistema d'equacions amb una matriu triangular superior, rebent una de triangular inferior passada per paràmetre (que transposada seria la triangular superior que volem), i un vector b passat com a paràmetre*/
int resLsup(int n, double **L, double *b, double *x, double tol){
	int i, k;
	double numerator, denominator, sum;
	numerator = b[n - 1];
	if (fabs(numerator) < tol) numerator = 0.;
	denominator = L[n - 1][n - 1];
	/* si el denominador és nul, el sistema és incompatible i no pot ser resolt*/
	if (fabs(denominator) < tol) return 1;
	x[n - 1] = numerator / denominator;
	for (i = n - 2; i >= 0; i--){
		sum = 0.;
		for (k = i + 1; k < n; k++){
			sum += L[k][i] * x[k];
		}
		numerator = b[i] - sum;
		if (fabs(numerator) < tol) numerator = 0.;
		denominator = L[i][i];
		/* si el denominador és nul, el sistema és incompatible i no pot ser resolt*/
		if (fabs(denominator) < tol) return 1;
		x[i] = numerator / denominator;
		if (fabs(x[i]) < tol) x[i] = 0.;
	}
	return 0;
}
