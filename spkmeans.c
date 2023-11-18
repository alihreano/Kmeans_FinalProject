#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "utils.h"
#include "spkmeans.h"

int main(int argc, char *argv[])
{
    FILE *f;
    int rows, cols;
    double **data;
    double **wam, **ddg, **gl;
    eigenHolder **jacobiRes;
    char c;
    char *goal = argv[1];
    char *inputFileName = argv[2];

    int i = 0, j, read;
    char *cord, *line = NULL;
    size_t len = 0;

    f = fopen(inputFileName, "r");
    if (argc == 0)
    {
        printf("An error has occurred\n");
        exit(EXIT_FAILURE);
    }
    
    /* calculating number of points */
    rows = 0;
    for (c = getc(f); c != EOF; c = getc(f))
    {
        if (c == '\n')
        {
            rows++;
        }
    }

    /* calculating the points dimension */
    rewind(f);
    cols = 1;
    for (c = getc(f); c != '\n'; c = getc(f))
    {
        if (c == ',')
        {
            cols++;
        }
    }
    rewind(f);

    /* reading the actual data to double** */
    data = (double **)malloc(rows * sizeof(double *));
    while ((read = getline(&line, &len, f)) != -1)
    {
        data[i] = (double *)malloc(cols * sizeof(double));
        if (data[i] == NULL)
        {
            printf("An error has occurred\n");
            exit(EXIT_FAILURE);
        }

        cord = strtok(line, ",");
        for (j = 0; j < cols; j++)
        {
            if (j == cols - 1)
            {
                cord = strtok(cord, "\n");
            }

            data[i][j] = strtod(cord, &cord + sizeof(double));
            cord = strtok(NULL, ",");
        }

        i++;
    }

    free(line);

    fclose(f);

    if (strcmp(goal, "wam") == 0)
    {
        wam = genWAM(data, rows, cols);
        printMatrix(wam, rows, rows);
        freeMatrix(wam, rows);
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        wam = genWAM(data, rows, cols);
        ddg = genDDG(wam, rows);
        printMatrix(ddg, rows, rows);
        freeMatrix(wam, rows);
        freeMatrix(ddg, rows);
    }
    else if (strcmp(goal, "gl") == 0)
    {
        wam = genWAM(data, rows, cols);
        ddg = genDDG(wam, rows);
        gl = genGL(wam, ddg, rows);
        printMatrix(gl, rows, rows);
        freeMatrix(wam, rows);
        freeMatrix(ddg, rows);
        freeMatrix(gl, rows);
    }
    else if (strcmp(goal, "jacobi") == 0)
    {
        jacobiRes = jacobi(data, rows);
        printEigs(jacobiRes, rows);
        freeEigenHolder(jacobiRes, rows);
    }
    else
    {
        printf("invalid goal");
    }

    freeMatrix(data, rows);

    return 0;
}