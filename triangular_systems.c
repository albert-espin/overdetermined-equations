/* Programa que calcula solucions de sistemes d'equacions sobredeterminats*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int resLinf(int, double **, double *, double *, double);
void prodMatVect(int, int, double **, double *, double *);
int resLsup(int, double **, double *, double *, double);

/* funció principal del programa*/
int main(void){
	int scanResult, n, i, j, veredict;
	double tol, remainderTerm, euclideanNorm;
	double * vector, * solutionVector, * productVector;
	double ** matrix;

	printf("Introdueix la tolerancia (per exemple 1e-3): ");
	scanResult = scanf("%le", &tol);

	/* obtenció de la dimensió*/
	printf("Introdueix la dimensio n: ");
	scanResult = scanf("%d", &n);

	/* obtenció de la matriu triangular inferior*/
	printf("Introdueix els elements de la part inferior (o diagonal) de la matriu (%dx%d)\n triangular inferior,separant-los amb espais:\n", n, n);
	matrix = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; i++){
		matrix[i] = (double*)malloc(n * sizeof(double));
		printf("Introdueix els %d elements de\n la part inferior (o diagonal) de la fila %d: ", i + 1, i + 1);
		/* registrem els valors de la part inferior o diagonal de la matriu*/
		for (j = 0; j <= i; j++){
			scanResult = scanf("%le", &matrix[i][j]);
		}
		/* omplim amb zeros la part superior de la matriu*/
		for (j = j; j < n; j++){
			matrix[i][j] = 0;
		}
	}

	/* mostra de la matriu*/
	if (matrix == NULL){
		printf("No existeix la matriu.\n");
		exit(1);
	}
	else{
		printf("La matriu correspon a:\n");
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				printf("%e ", matrix[i][j]);
			}
			printf("\n");
		}
	}


	/* obtenció del vector*/
	printf("Introdueix les %d components del vector: ", n);
	vector = (double*)malloc(n * sizeof(double));
	for (i = 0; i < n; i++){
		scanResult = scanf("%le", &vector[i]);
	}

	solutionVector = (double*)malloc(n * sizeof(double));
	veredict = resLinf(n, matrix, vector, solutionVector, tol);

	/*si el sistema no té solució ho indiquem*/
	if (veredict == 1){
		printf("El sistema Lx = b no te solucio\n");
	}
	/* si té la solució, prosseguim i la mostrem*/
	else{
		/* mostra de la solució del sistema entre la matriu i el vector donats*/
		printf("\nEl vector solucio resulta ser:\n");
		for (i = 0; i < n; i++){
			printf("%e\n", solutionVector[i]);
		}

		/* mostra del producte entre la matriu i el vector solució del sistema*/
		printf("\nAddicionalment, multipliquem la matriu inicial\npel vector de solucions del sistema.\n");
		productVector = (double*)malloc(n * sizeof(double));
		prodMatVect(n, n, matrix, solutionVector, productVector);
		printf("Aixo ens permet trobar la bondat del sistema, fent la norma euclidiana \ndels elements de ||L*x - b||\n");
		euclideanNorm = 0.;
		for (i = 0; i < n; i++){
			remainderTerm = fabs(productVector[i] - vector[i]);
			euclideanNorm += remainderTerm * remainderTerm;
		}
		euclideanNorm = sqrt(euclideanNorm);
		printf("%e\n", euclideanNorm);


		/* resolució del sistema d'equacions entre la transposta de la matriu inicial (per tant triangular superior) i el vector original b*/
		printf("\nAra transposem la matriu original i resolem el sistema\nentre aquesta i el vector inicial.\n");
		veredict = resLsup(n, matrix, vector, solutionVector, tol);

		/*si el sistema no té solució ho indiquem*/
		if (veredict == 1){
			printf("El sistema no te solucio\n");
		}

		else{
			printf("El vector solucio resulta ser:\n");
			for (i = 0; i < n; i++){
				printf("%e\n", solutionVector[i]);
			}
		}
	}
  
  return scanResult;
}
