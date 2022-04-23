#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

static inline int mindex(const int row, const int col, const int size);

void allocateMatrix(double **a, double **l, double **u, double **ver, double **aorigin, size_t NBYTES);
void generateMatrix(double *a, double *aorigin, int size);
void factorization(double *a, double *l, double *u, int size);
int verification(double *l, double *u, double *aorigin, double *ver, int size);
void printMatrix(double *matrix, int size);


clock_t start_point, end_point;

#define FILENAME "output.txt"

int main(int argc, char * argv[]){

	if(argc < 3) {
		printf("<size> <threadNumber>!\n");
		return EXIT_FAILURE;
	}
	
	
	const int size  = atoi(argv[1]);
	const size_t NBYTES = size * size * sizeof(double);

	const int threadNumber = atoi(argv[2]);

	omp_set_num_threads(threadNumber);

	double * a, *l, *u, *ver, *aorigin;

	allocateMatrix(&a, &l, &u, &ver, &aorigin, NBYTES);
	generateMatrix(a,aorigin, size);

	factorization(a, l, u, size);


	
	FILE * fp = fopen(FILENAME, "a+");

	float operating_time = (float) (end_point - start_point) / CLOCKS_PER_SEC ;

	int isVerified =  verification(l, u, aorigin, ver, size);

	if(isVerified) {
		fprintf(fp, "%d ", size);
		fprintf(fp, "%f ", operating_time);
		fprintf(fp, "%s\n", "VERIFIED");
	}

	else {
		fprintf(fp, "%d ", size);
		fprintf(fp, "%f ", operating_time);
		fprintf(fp, "%s\n", "NOT VERIFIED");
	}

	fclose(fp);


	free(a);
	free(l);
	free(u);
	free(ver);
	free(aorigin);

	return EXIT_SUCCESS;

}


void printMatrix(double *matrix, int size) {
	unsigned int i, j;

	for (i = 0; i < size; i++) {
		for(j = 0; j < size; j++){

				printf("| %f | ", matrix[mindex(i, j, size)]);
		}
		printf("\n");
	}
}

void allocateMatrix(double **a, double **l, double **u, double **ver, double **aorigin, size_t NBYTES) {
	*a = malloc(NBYTES);	
	*l = malloc(NBYTES);
	*u = malloc(NBYTES);
	*ver = malloc(NBYTES);
	*aorigin = malloc(NBYTES);
}


void generateMatrix(double *a, double *aorigin, int size) {
	
	for(unsigned int i = 0;  i < size; i++){
		for(unsigned int j = 0; j < size; j++){
			 double val =  (rand() % 5 + 1);
			 a[mindex(i,j, size)] = val;
			 aorigin[mindex(i,j, size)] = val;
		}
	}
}



void factorization(double *a, double *l, double *u, int size) {
	start_point = clock();

	unsigned int i, k, j;

	#pragma omp parallel shared(a, l, u) private(i, j, k)
	{
		for(k=0; k<size; k++){
        	for(i=k+1; i<size; i++) {
 				//a[i][k] /= a[k][k];
        		a[mindex(i, k, size)] /= a[mindex(k, k, size)];	
        	}
        	

        	#pragma omp barrier
        	#pragma omp for schedule(static, 1)
        	for(i=k+1; i<size; i++){
            	for(j=k+1; j<size; j++){
               		//a[i][j] -= a[i][k] * a[k][j];
                	a[mindex(i, j, size)] -= a[mindex(i, k, size)] * a[mindex(k, j, size)]; 
            	}
        	}   
        }	

        #pragma omp barrier
        #pragma omp for schedule(static, 1)
		for(i=0; i<size; i++){
        	for(j=0; j<=i; j++){
            	if(i==j) {
            		//l[i][j] = 1.0f;
            		l[mindex(i, j, size)] = 1.0; 
            	}   
            	else  {
            		//l[i][j] = a[i][j];
            		l[mindex(i, j, size)] =  a[mindex(i, j, size)]; 	
            	}      
        	}
    	}

    	#pragma omp for schedule(static, 1)
    	for(i=0; i < size; i++){
        	for(j=i; j < size; j++){
            	//u[i][j] = a[i][j];
            	u[mindex(i, j, size)] = a[mindex(i, j, size)];
        	}
    	}

	}
    

	end_point = clock();
}

int verification(double *l, double *u, double *aorigin, double *ver, int size) {

	int isVerified = 1;

	unsigned int i, j, k;
	for (i = 0; i < size; i++) {
		for(j = 0; j < size; j++ ) {
			//ver[i][j] = 0.0f;
			ver[mindex(i, j, size)] = 0.0;
			for(k = 0; k < size; k++) {
				//ver[i][j] += l[i][k] * u[k][j];
				ver[mindex(i,j, size)] += l[mindex(i, k, size)] * u[mindex(k, j, size)]; 
			}

		}
	}

	for (i=0; i < size; i++) {
		for(j = 0; j < size; j++){
			if(fabs(ver[mindex(i, j, size)] - aorigin[mindex(i, j, size)]) > 0.001){
				isVerified = 0;
				break;
			}	

		}
	}

	return isVerified;

}


static inline int mindex(const int row, const int col, const int size) {
	return row * size + col;	
}
