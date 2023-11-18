
#include "kmeanspp.h"


void print_res(double **curr, int k, int d)
{
    int i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (j == (d - 1))
            {
                printf("%.4f", curr[i][j]);
                printf("\n");
            }
            else
            {
                printf("%.4f", curr[i][j]);
                printf("%s", ",");
            }
        }
    }
}

double **empty_matrix(int rows, int cols)
{
    int i, j;
    double **matrix = (double **)calloc(rows, sizeof(double *));
    if (matrix == NULL)
    {
        return NULL;
    }
    for (i = 0; i < rows; i++)
    {
        matrix[i] = (double *)calloc(cols, sizeof(double));
        if (matrix[i] == NULL)
        {
            for (j = 0; j < i; j++)
            {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

double **init_Kclusters(double **mat, int k, int d)
{
    double **res;
    int i, j;
    res = empty_matrix(k, d);
    if (res == NULL)
    {
        return NULL;
    }
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            res[i][j] = mat[i][j];
        }
    }
    return res;
}

double **update_centroids(double ***groups, int *points_num, int k, int d)
{
    double **res;
    int i, j, l;
    double avg, x;
    res = (double **)calloc(k, sizeof(double *));
    if (res == NULL)
    {
        printf("An Error Has Occurred\n");
        return NULL;
    }
    x = 0.0;
    for (i = 0; i < k; i++)
    {
        res[i] = (double *)calloc(d, sizeof(double));
        if (res[i] == NULL)
        {
            printf("An Error Has Occurred\n");
            return NULL;
        }
        for (j = 0; j < d; j++)
        {
            x = 0.0;
            for (l = 0; l < points_num[i]; l++)
            {
                x += groups[i][l][j];
            }
            avg = x / points_num[i];
            res[i][j] = avg;
        }
    }
    return res;
}

double ***make_empty_groups(int k, int *points_num, int d)
{
    int i, j, n;
    double ***res = (double ***)calloc(k, sizeof(double **));
    if (res == NULL)
    {
        return NULL;
    }
    for (i = 0; i < k; i++)
    {
        n = points_num[i];
        res[i] = (double **)calloc(n, sizeof(double *));
        if (res[i] == NULL)
        {
            return NULL;
        }
        for (j = 0; j < n; j++)
        {
            res[i][j] = (double *)calloc(d, sizeof(double));
            if (res[i][j] == NULL)
            {
                return NULL;
            }
        }
    }
    return res;
}
double Euclidean_Distance(double *vec1, double *vec2, int d)
{
    double res = 0.0;
    int i;
    for (i = 0; i < d; i++)
    {
        res = res + (double)pow((double)(vec1[i]) - (double)(vec2[i]), 2.0);
    }
    return (double)sqrt(res);
}
bool check_convergence(int epsilon, double **curr, double **prev, bool first_iter, int k)
{
    int i;
    if (first_iter)
    {
        return false;
    }
    else
    {
        for (i = 0; i < k; i++)
        {
            if (Euclidean_Distance(curr[i], prev[i], d) >= epsilon)
            {
                return false;
            }
        }
        return true;
    }
    return false;
}
/*
void freeMemory(double** matrix ,int len){
        int i;
        if(matrix == NULL){
            return;
        }
        for(i = 0; i < len ; i++){
            if(matrix[i] == NULL){
                continue;
            }
            free(matrix[i]);
        }
        free(matrix);
}
*/

void freeMemory(double **matrix, int len)
{
    int i;
    for (i = 0; i < len; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void free_empty_groups(double ***groups, int k, int *points_num)
{
    int i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < points_num[i]; j++)
        {
            free(groups[i][j]);
        }
        free(groups[i]);
    }
    free(groups);
}

double **kmeans(double **mat, double **curr_clusters, int n, int k, int d, int iter, double eps)
{

    /*int i;*/

    /* int count=0 ;*/
    double **prev_clusters;
    bool first_iter;
    int iter_count;
    double ***groups = NULL;
    int x_i;
    int min_index = 0;
    int cluster_ind;
    double min_dist;
    double curr_dist;
    int x;
    int *points_num;
    int *ind_vec;
    if ((k <= 1) || (k >= n) || (k != floor(k)))
    {
        printf("Invalid number of clusters!");
        printf("\n");
        return NULL;
    }
    if ((iter <= 1) || (iter >= 1000) || (iter != floor(iter)))
    {
        printf("Invalid maximum iteration!");
        printf("\n");
        return NULL;
    }

    prev_clusters = empty_matrix(k, d);

    if (prev_clusters == NULL)
    {
        printf("An Error Has Occurred");
        printf("\n");
        return NULL;
    }
    first_iter = true;
    iter_count = 0;
    points_num = (int *)calloc(k, sizeof(int));
    if (points_num == NULL)
    {
        printf("An Error Has Occurred");
        printf("\n");
        return NULL;
    }
    ind_vec = (int *)calloc(k, sizeof(int));
    if (ind_vec == NULL)
    {
        printf("An Error Has Occurred");
        printf("\n");
        return NULL;
    }
    x = 0;
    while ((!check_convergence(eps, curr_clusters, prev_clusters, first_iter, k)) && (iter_count < iter))
    {
        first_iter = false;
        prev_clusters = curr_clusters;
        for (x_i = 0; x_i < n; x_i++)
        {
            min_dist = -1;
            for (cluster_ind = 0; cluster_ind < k; cluster_ind++)
            {
                curr_dist = Euclidean_Distance(mat[x_i], curr_clusters[cluster_ind], d);
                if ((curr_dist < min_dist) || (min_dist == -1))
                {
                    min_dist = curr_dist;
                    min_index = cluster_ind;
                }
            }
            points_num[min_index]++;
            ind_vec[min_index]++;
        }
        groups = make_empty_groups(k, points_num, d);
        if (groups == NULL)
        {
            printf("An Error Has Occurred");
            printf("\n");
            return NULL;
        }
        for (x_i = 0; x_i < n; x_i++)
        {
            min_dist = -1;
            for (cluster_ind = 0; cluster_ind < k; cluster_ind++)
            {
                curr_dist = Euclidean_Distance(mat[x_i], curr_clusters[cluster_ind], d);
                if ((curr_dist < min_dist) || (min_dist == -1))
                {
                    min_dist = curr_dist;
                    min_index = cluster_ind;
                }
            }
            x = ind_vec[min_index] - 1;
            groups[min_index][x] = mat[x_i];
            ind_vec[min_index] = ind_vec[min_index] - 1;
        }
        curr_clusters = update_centroids(groups, points_num, k, d);
        points_num = (int *)calloc(k, sizeof(int));
        if (points_num == NULL)
        {
            printf("An Error Has Occurred");
            printf("\n");
            return NULL;
        }
        iter_count++;
    }

    free_empty_groups(groups, k, points_num);
    free(points_num);
    freeMemory(prev_clusters, k);
    return curr_clusters;
}
