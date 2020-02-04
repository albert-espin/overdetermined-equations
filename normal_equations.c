/* Programa que calcula solucions de sistemes d'equacions sobredeterminats*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int resLinf(int, double **, double *, double *, double);
int resLsup(int, double **, double *, double *, double);
int ldlt(int, double **, double *, double);
void prodMatVect(int, int, double **, double *, double *);
void prodMatMat(int, int, int, double **, double **, double **);

/* funció principal del programa*/
int main(void){
	int scanResult, n, m, i, j, decompositionVeredict;
	double tol, remainderTerm, euclideanNorm;
	double * vector, *solutionVector, * productVector;
	double ** matrix, ** transpose, ** productMatrix;

	printf("Trobarem una solucio de les infinites possibles a un sistema sobredeterminat.\n");

	printf("Introdueix la tolerancia (per exemple 1e-3): ");
	scanResult = scanf("%le", &tol);

	/* obtenció de le dimensions*/
	printf("Introdueix la dimensio m (nombre de files): ");
	scanResult = scanf("%d", &m);
	printf("Introdueix la dimensio n (nombre de columnes): ");
	scanResult = scanf("%d", &n);

	/* els sistemes sobredeterminats exigeixen matrius amb més files que columnes; si això no es compleix no estem davant d'un sistem sobredeterminat*/
	if (n >= m){
		printf("Amb aquestes dades no podrem tenir un sistema sobredeterminat.\nFinalitzem l'execucio.");
		exit(0);
	}

	/* obtenció de la matriu*/
	printf("Ara volem obtenir els (%dx%d) elements\nde la matriu, separats amb espais:\n", m, n);
	matrix = (double**)malloc(m * sizeof(double*));
	for (i = 0; i < m; i++){
		matrix[i] = (double*)malloc(n * sizeof(double));
		printf("Introdueix els %d elements de la fila %d: ", n, i + 1);
		/* registrem els valors de cada fila de la matriu*/
		for (j = 0; j < n; j++){
			scanResult = scanf("%le", &matrix[i][j]);
		}
	}

	/* obtenció de la transposta*/
	transpose = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; i++){
		transpose[i] = (double*)malloc(m * sizeof(double));
		for (j = 0; j < m; j++){
			transpose[i][j] = matrix[j][i];
		}
	}

	/* mostra de les matrius*/
	if (matrix == NULL || transpose == NULL){
		printf("No existeix la matriu.\n");
		exit(1);
	}
	else{
		printf("La matriu correspon a:\n");
		for (i = 0; i < m; i++){
			for (j = 0; j < n; j++){
				printf("%e ", matrix[i][j]);
			}
			printf("\n");
		}

		printf("La transposta correspon a:\n");
		for (i = 0; i < n; i++){
			for (j = 0; j < m; j++){
				printf("%e ", transpose[i][j]);
			}
			printf("\n");
		}
	}

	/* obtenció del producte At*A*/
	productMatrix = (double**)malloc(m * sizeof(double*));
	prodMatMat(n, m, n, transpose, matrix, productMatrix);
	printf("\nEl producte de la transposta per\nla matriu inicial dona la matriu:\n");
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			printf("%e ", productMatrix[i][j]);
		}
		printf("\n");
	}

	/* obtenció del vector b*/
	printf("\nAra volem un vector de termes indepentens.\nIntrodueix les %d components del vector: ", m);
	vector = (double*)malloc(m * sizeof(double));
	for (i = 0; i < m; i++){
		scanResult = scanf("%le", &vector[i]);
	}

	/* obtenció del producte At*b*/
	solutionVector = (double*)malloc(m * sizeof(double));
	prodMatVect(n, m, transpose, vector, solutionVector);
	printf("\nTrobem ara el vector del producte de la transposta pel vector.\n");

	/* descomposició LDLt completa */
	decompositionVeredict = ldlt(n, productMatrix, solutionVector, tol);
	if (decompositionVeredict == 0){
		/* mostra de la solució */
		printf("\nLa solucio del sistema Ax* = b ((At*A)x) = At*b)\namb descomposicio LDLt correspon a:\n");
		for (i = 0; i < n; i++){
			printf("%e \n", solutionVector[i]);
		}

		/* mostra del producte entre la matriu i el vector solució del sistema*/
		printf("\nAddicionalment, multipliquem la matriu inicial\npel vector de solucions del sistema.\n");
		productVector = (double*)malloc(n * sizeof(double));
		prodMatVect(n, n, matrix, solutionVector, productVector);
		printf("Aixo ens permet trobar la bondat del sistema, fent la norma euclidiana \ndels elements de ||A*x - b||:\n");
		euclideanNorm = 0.;
		for (i = 0; i < n; i++){
			remainderTerm = fabs(productVector[i] - vector[i]);
			euclideanNorm += remainderTerm * remainderTerm;
		}
		euclideanNorm = sqrt(euclideanNorm);
		printf("%e\n", euclideanNorm);
	}

	else{
		printf("\nEl sistema Ax = b no te solucio.");
	}

	return scanResult;
}

/*funció que donades unes matrius A i B de dimensions respectives m x n i n x p calcula una matriu producte C de dimensions n x p*/
void prodMatMat(int m, int n, int p, double **A, double **B, double **C){
	int i, j, k;
	double sum;
	for (i = 0; i < m; i++){
		C[i] = (double *)malloc(p * sizeof(double));
		for (j = 0; j < p; j++){
			sum = 0;
			for (k = 0; k < n; k++){
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
		}
	}
}



/* funció que donada una matriu simètrica A de dimensions n x n i un vector de termes independents v troba la solució al sistema A*x = v descomposant primer A com a L*D*Ltransposada, tenint en compte una tolerència tol i retornant 0 si el sistema té solució determinada i 1 en cas contrari*/
int ldlt(int n, double **A, double *v, double tol){
	int i, j, k, veredict;
	double sum;
	double *y, *z;
	double **L, **D;

	/* reserva de memòria per a D i L; omplim D de zeros a tot arreu menys a la diagonal i fem 1 els termes de la diagonal de L, mentre que els termes per sobre de la diagonal seran 0*/
	D = (double**)malloc(n * sizeof(double*));
	L = (double**)malloc(n * sizeof(double*));
	for (i = 0; i < n; i++){
		D[i] = (double*)malloc(n * sizeof(double));
		for (j = 0; j < n; j++){
			if (i != j) D[i][j] = 0;
		}
		L[i] = (double*)malloc(n * sizeof(double));
		L[i][i] = 1;
		for (j = i + 1; j < n; j++){
			L[i][j] = 0;
		}
	}

	/* obtenció de D i L */
	for (k = 0; k < n; k++){

		/* càlcul de dkk */
		sum = 0;
		for (j = 0; j < k; j++){
			sum += L[k][j] * L[k][j] * D[j][j];
		}
		D[k][k] = A[k][k] - sum;
		/* el sistema és indeterminat si hi ha una fila de zeros*/
		if (fabs(D[k][k]) < tol) return 1;

		/* càlcul de lik */
		for (i = k + 1; i < n; i++){
			sum = 0;
			for (j = 0; j < k; j++){
				sum += L[i][j] * L[k][j] * D[j][j];
			}
			L[i][k] = (A[i][k] - sum) / D[k][k];
			if (fabs(L[i][k]) < tol) L[i][k] = 0;
		}
	}

	printf("\nLa matriu D correspon a:\n");
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			printf("%e ", D[i][j]);
		}
		printf("\n");
	}

	printf("\nLa matriu L correspon a:\n");
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			printf("%e ", L[i][j]);
		}
		printf("\n");
	}

	/* L * z = b */
	z = (double*)malloc(n * sizeof(double));
	veredict = resLinf(n, L, v, z, tol);
	if (veredict == 1) return 1;

	/* D * y = z */
	y = (double*)malloc(n * sizeof(double));
	veredict = resLinf(n, D, z, y, tol);
	if (veredict == 1) return 1;

	/* Lt * x = y */
	veredict = resLsup(n, L, y, v, tol);
	if (veredict == 1) return 1;

	return 0;
}
