#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "math.h"
#include <float.h>
#include <stdbool.h>
#include <math.h>
#define assert__(x) for (; !(x); assert(x))
#define PY_SSIZE_T_CLEAN

int iter;
int d;
int n;
int k;
double **empty_matrix(int rows, int cols);
double **init_Kclusters(double **mat, int k, int d);
void freeMemory(double **matrix, int len);
bool check_convergence(int epsilon, double **curr, double **prev, bool first_iter, int k);
double Euclidean_Distance(double *vec1, double *vec2, int d);
double ***make_empty_groups(int k, int *points_num, int d);
double **update_centroids(double ***groups, int *points_num, int k, int d);
void print_res(double **curr, int k, int d);
void free_empty_groups(double ***groups, int k, int *points_num);
double **kmeans(double **mat, double **curr_clusters, int n, int k, int d, int iter, double eps);