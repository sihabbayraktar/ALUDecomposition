#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

void create(double ***a, double ***l, double ***u, double ***ver, double ***atemp, int size);
void decompose(double **a, double **l, double **u, int size);
void fill(double **a, double **u, double **atemp, int size);
int correctLU(double **l, double **u, double **atemp, double **ver, int size);
void destroy(double **matrix, int size); 

#define TOLERABLE_ERROR 1e-03

int main(int argc, char * argv[]){

	if(argc < 2) {
		printf("./sequential <size>!\n");
		return 1;
	}
	double **a, **l, **u, **ver, **atemp; 
	struct timeval start_point, end_point;
    double elapsed_time;

	int size  = atoi(argv[1]);
	

	create(&a, &l, &u, &ver, &atemp, size);
	fill(a, u, atemp, size);


	gettimeofday(&start_point, NULL);
	decompose(a, l, u, size);
	gettimeofday(&end_point, NULL);
	
	elapsed_time = (end_point.tv_sec - start_point.tv_sec) * 1e6;
    elapsed_time = (elapsed_time + (end_point.tv_usec -  start_point.tv_usec)) * 1e-6;

	int correct = correctLU(l, u, atemp, ver, size);


	
	FILE * fp = fopen("result.txt", "a+");

	if(correct) {
		fprintf(fp, "%d ", size);
		fprintf(fp, "%.6f ", elapsed_time);
		fprintf(fp, "%s\n", "MATCH");
	}

	else {
		fprintf(fp, "%d ", size);
		fprintf(fp, "%.6f ", elapsed_time);
		fprintf(fp, "%s\n", "NOT MATCH");
	}

	fclose(fp);


	destroy(a, size);
	destroy(l, size);
	destroy(u, size);
	destroy(ver, size);
	destroy(atemp, size);

	return 0;
}


void destroy(double **matrix, int size) {
	for(int i = 0; i < size; i++){
		free(matrix[i]);
	}
	free(matrix);
}


void create(double ***a, double ***l, double ***u, double ***ver, double ***atemp, int size) {

	*a = (double **)malloc(sizeof(double *) * size);
	*l = (double **)malloc(sizeof(double *) * size);
	*u = (double **)malloc(sizeof(double *) * size);
	*ver = (double **)malloc(sizeof(double *) * size);
	*atemp =  (double **)malloc(sizeof(double *) * size);

	for(int i = 0; i < size; i++){
		(*a)[i] = (double *) malloc(sizeof(double) * size);
		(*l)[i] = (double *) malloc(sizeof(double) * size);
		(*u)[i] = (double *) malloc(sizeof(double) * size);
		(*ver)[i] = (double *) malloc(sizeof(double) * size);
		(*atemp)[i] = (double *) malloc(sizeof(double) * size);


		for(int j = 0; j < size; j++){
			(*a)[i][j] = 0.0f;
			(*l)[i][j] = 0.0f;
			(*u)[i][j] = 0.0f;
			(*ver)[i][j] = 0.0f;
			(*atemp)[i][j] = 0.0f;
		}
	}	
}

void fill(double **a, double **u, double **atemp, int size) {
	
	for(int i = 0;  i < size; i++){
		for(int j = 0; j < size; j++){
			 double val = (double) (rand() % 10 + 1);
			 a[i][j] = val;
			 atemp[i][j] = val;
		}
	}
}


int correctLU(double **l, double **u, double **atemp, double **ver, int size) {

	int correct = 1;

	
	for (int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++ ) {
			ver[i][j] = 0.0;
			for(int k = 0; k < size; k++) {
				ver[i][j] += l[i][k] * u[k][j]; 
			}

		}
	}

	for (int i=0; i < size; i++) {
		for(int j = 0; j < size; j++){
			if(fabs(atemp[i][j] -ver[i][j]) > TOLERABLE_ERROR){
				correct = 0;
				break;
			}	

		}
	}

	return correct;

}

void decompose(double **a, double **l, double **u, int size) {

	int index, index2, index3;

    for(index=0; index<size; index++){
        for(index2=index+1; index2<size; index2++){
            a[index2][index] /= a[index][index];
        }
        for(index2=index+1; index2<size; index2++){
            for(index3=index+1; index3<size; index3++){
                a[index2][index3] -= a[index2][index] * a[index][index3];
            }
        }
    }

    for(index2=0; index2<size; index2++){
        for(index3=0; index3<=index2; index3++){
            if(index2==index3)    l[index2][index3] = 1.0f;
            else        l[index2][index3] = a[index2][index3];
        }
    }

    for(index2=0; index2 < size; index2++){
        for(index3=index2; index3 < size; index3++){
            u[index2][index3] = a[index2][index3];
        }
    }
}


