#include "utils.h"

/**
 * calculates the Weighted Adjacency Matrix of a given matrix.
 * by the formula W_ij = exp(-(||x_i - x_j||^2) / 2) for all i != j
 * and W_ii = 0
 */
double **genWAM(double **matrix, int rows, int cols)
{
    double **weightedAdjMatrix = (double **)calloc(rows, sizeof(double *));
    int i, j;

    for (i = 0; i < rows; i++)
    {
        weightedAdjMatrix[i] = (double *)calloc(rows, sizeof(double));
        for (j = 0; j < rows; j++)
        {
            if (j != i)
            {
                weightedAdjMatrix[i][j] = vectorsWeight(matrix[i], matrix[j], cols);
            }
            else
            {
                weightedAdjMatrix[i][j] = 0;
            }
        }
    }

    return weightedAdjMatrix;
}

/**
 * calculates the Diagonal Degree Matrix from given weighted adjacency matrix
 * returns diagonal matrix with sum(W_ij) [j=1, j=n] in the entry ii
 */
double **genDDG(double **wam, int n)
{
    double **ddg = calloc(n, sizeof(double));
    int i, j, k;
    double sum;

    for (i = 0; i < n; i++)
    {
        ddg[i] = calloc(n, sizeof(double));
        for (j = 0; j < n; j++)
        {
            if (j == i)
            {
                sum = 0;
                for (k = 0; k < n; k++)
                {
                    sum += wam[i][k];
                }
                ddg[i][j] = sum;
            }
            else
            {
                ddg[i][j] = 0;
            }
        }
    }

    return ddg;
}

/**
 * simply subtract matrices,
 * used to calculate the grpah laplacian by subtracting the WAM from DDG
 */
double **genGL(double **wam, double **ddg, int n)
{
    double **gl = calloc(n, sizeof(double));
    int i, j;

    for (i = 0; i < n; i++)
    {
        gl[i] = calloc(n, sizeof(double));
        for (j = 0; j < n; j++)
        {
            gl[i][j] = ddg[i][j] - wam[i][j];
        }
    }

    return gl;
}
/**
 * return exp(-(||u - v||^2) / 2) for given u, v double arrays
 */
double vectorsWeight(double *v, double *u, int cols)
{
    double dist = 0;
    int i;

    for (i = 0; i < cols; i++)
    {
        dist += pow((v[i] - u[i]), 2.0);
    }

    return exp(-1 * dist / 2);
}

/* -----jacobi----- */

/**
 * find what is the max(abs(off-diagonal elements)) of given matrix,
 * updates the result into given struct
 */
/* TODO: change the struct? */
void getMaxOffDiagonal(double **matrix, int rows, rotationMatrixData *rotationMatrix)
{
    int i, j;
    rotationMatrix->i = 0;
    rotationMatrix->j = 1;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < rows; j++)
        {
            /* we looking for the off diagonal largest number */
            if (i == j)
            {
                continue;
            }

            if (fabs(matrix[i][j]) > fabs(matrix[rotationMatrix->i][rotationMatrix->j]))
            {
                rotationMatrix->i = i;
                rotationMatrix->j = j;
            }
        }
    }
}

/**
 * calculates phi according to: θ = (A_jj - A_ii) / 2A_ij
 */
double calcPhi(double **matrix, rotationMatrixData *rotationMatrix)
{
    return (matrix[rotationMatrix->j][rotationMatrix->j] - matrix[rotationMatrix->i][rotationMatrix->i]) / (2 * matrix[rotationMatrix->i][rotationMatrix->j]);
}

/**
 * calculates T according to: t = sign(θ) / (|θ| +√(θ^2 + 1)) from given phi
 */
double calcT(double phi)
{
    return sign(phi) / (fabs(phi) + sqrt(pow(phi, 2.0) + 1));
}

/**
 * calculates c according to c = 1 / √(t^2 + 1) from given t
 */
double calcC(double t)
{
    return 1 / sqrt(pow(t, 2.0) + 1);
}

/**
 * returns sign of double value.
 * Note: we defined sign(0) to be 1.
 */
int sign(double val)
{
    return val >= 0 ? 1 : -1;
}

/**
 * returns an identity squared matrix of shape = (rows, rows)
 */
double **getIdMatrix(int rows)
{
    double **id;
    int i;

    id = malloc(rows * sizeof(double *));
    if (id == NULL)
    {
        printf("An error has occurred\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < rows; i++)
    {
        id[i] = calloc(rows, sizeof(double));
        id[i][i] = 1.0;
    }

    return id;
}

/**
 * returns the matrix AB from given matrices A, B
 */
double **matrixProduct(double **matrixA, double **matrixB, int rows)
{
    int i, j, k;
    double **result = malloc(rows * sizeof(double *));
    double resAcc;
    for (i = 0; i < rows; i++)
    {
        result[i] = malloc(rows * sizeof(double));

        for (j = 0; j < rows; j++)
        {
            resAcc = 0;
            for (k = 0; k < rows; k++)
            {
                resAcc += matrixA[i][k] * matrixB[k][j];
            }
            result[i][j] = resAcc;
        }
    }

    return result;
}

/**
 * calculate the uptadetd A' according to the insturctions in section 6
 */
double **updateAtag(double **matrix, rotationMatrixData *rotationMatrix, int rows)
{
    int i, j;
    double **updatedMatrix;

    updatedMatrix = (double **)malloc(rows * sizeof(double *));
    if (updatedMatrix == NULL)
    {
        printf("An error has occurred\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < rows; i++)
    {
        updatedMatrix[i] = malloc(rows * sizeof(double));
        if (updatedMatrix[i] == NULL)
        {
            printf("An error has occurred\n");
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < rows; j++)
        {
            if (i == rotationMatrix->i && j == rotationMatrix->i)
            { /* a'_ii = c^2 * a_ii + s^2 * a_jj − 2 * s * c * a_ij */

                updatedMatrix[i][j] = pow(rotationMatrix->c, 2) * matrix[rotationMatrix->i][rotationMatrix->i] +
                                      pow(rotationMatrix->s, 2) * matrix[rotationMatrix->j][rotationMatrix->j] - 2 * rotationMatrix->s * rotationMatrix->c * matrix[rotationMatrix->i][rotationMatrix->j];
            }
            else if (i == rotationMatrix->j && j == rotationMatrix->j)
            { /* a'_jj = s^2 * a_ii + c^2 * a_jj − 2 * s * c * a_ij */
                updatedMatrix[i][j] = pow(rotationMatrix->c, 2) * matrix[rotationMatrix->j][rotationMatrix->j] + pow(rotationMatrix->s, 2) * matrix[rotationMatrix->i][rotationMatrix->i] +
                                      2 * rotationMatrix->s * rotationMatrix->c * matrix[rotationMatrix->i][rotationMatrix->j];
            }
            else if ((i == rotationMatrix->j && j == rotationMatrix->i) || (i == rotationMatrix->i && j == rotationMatrix->j))
            { /* a'_jj = (c^2 - s^2) * a_ij + s * c *(a_ii - a_jj) => a'_ij = 0 */
                updatedMatrix[i][j] = 0;
            }
            else if (j == rotationMatrix->i)
            { /* a'_ir */
                updatedMatrix[i][j] = rotationMatrix->c * matrix[i][rotationMatrix->i] - rotationMatrix->s * matrix[i][rotationMatrix->j];
            }
            else if (i == rotationMatrix->i)
            { /* a'_ri = c * a_ri - s * a_rj      r != i, j */

                updatedMatrix[i][j] = rotationMatrix->c * matrix[rotationMatrix->i][j] - rotationMatrix->s * matrix[rotationMatrix->j][j];
            }
            else if (j == rotationMatrix->j)
            { /* a'_rj = = c * a_rj - s * a_ri      r != i, j */

                updatedMatrix[i][j] = rotationMatrix->c * matrix[i][rotationMatrix->j] + rotationMatrix->s * matrix[i][rotationMatrix->i];
            }
            else if (i == rotationMatrix->j)
            { /* a'_jr */
                updatedMatrix[i][j] = rotationMatrix->c * matrix[rotationMatrix->j][j] + rotationMatrix->s * matrix[rotationMatrix->i][j];
            }
            else
            {
                updatedMatrix[i][j] = matrix[i][j];
            }
        }
    }

    return updatedMatrix;
}

/**
 * calculate rotation matrix P according to the instructions in section 2
 */
double **genRotationMatrix(rotationMatrixData *rotationMatrix, int rows)
{
    int i, j;
    double **res = malloc(rows * sizeof(double *));
    if (res == NULL)
    {
        printf("An error has occurred\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < rows; i++)
    {
        res[i] = malloc(rows * sizeof(double));
        for (j = 0; j < rows; j++)
        {
            if (i == j && (i == rotationMatrix->i || j == rotationMatrix->j))
            {
                res[i][j] = rotationMatrix->c;
                continue;
            }

            if (i == j)
            {
                res[i][j] = 1;
                continue;
            }

            if ((i == rotationMatrix->i && j == rotationMatrix->j))
            {
                res[i][j] = rotationMatrix->s;
                continue;
            }

            if ((j == rotationMatrix->i && i == rotationMatrix->j))
            {
                res[i][j] = -1 * rotationMatrix->s;
                continue;
            }

            res[i][j] = 0;
        }
    }

    return res;
}

/*
 * calculate the *squared* 2-norm of the elements off the diagonal
 */
double offDiagEuclidDistance(double **matrix, int rows)
{
    double sumSquare = 0;
    int i, j;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < rows; j++)
        {
            if (i != j)
            {
                sumSquare += pow(matrix[i][j], 2.0);
            }
        }
    }

    return sumSquare;
}

/**
 * takes the eigenvalues from diagonal matrix and eigvectors from another matrix,
 * returns them as eigenHolder struct
 */
eigenHolder **matrixToEig(double **matrix, double **eigenVectors, int rows)
{
    int i, j;
    double *eigenVector;
    eigenHolder **eigs = malloc(rows * sizeof(eigenHolder *));
    if (eigs == NULL)
    {
        printf("An error has occurred\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < rows; i++)
    {
        eigenHolder *eigenParams = malloc(sizeof(eigenHolder));
        if (eigenParams == NULL)
        {
            printf("An error has occurred\n");
            exit(EXIT_FAILURE);
        }

        eigenParams->eigValue = matrix[i][i];
        eigenVector = malloc(rows * sizeof(double));
        if (eigenVector == NULL)
        {
            printf("An error has occurred\n");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < rows; j++)
        {
            eigenVector[j] = eigenVectors[j][i];
        }
        eigenParams->eigenVector = eigenVector;
        eigs[i] = eigenParams;
    }

    return eigs;
}

/**
 * performs jacobi algorithm as described in the assignment
 */
eigenHolder **jacobi(double **matrix, int rows)
{
    double **updatedMatrix, **updatedEigVectors, **rotMatrix;
    int isConverged = 0, iterCounter = 0;
    eigenHolder **eigs;
    rotationMatrixData *rotationMatrix = malloc(sizeof(rotationMatrixData));
    double **eigenVectors = getIdMatrix(rows);

    while (!isConverged)
    {
        iterCounter++;
        getMaxOffDiagonal(matrix, rows, rotationMatrix);

        rotationMatrix->phi = calcPhi(matrix, rotationMatrix);
        rotationMatrix->t = calcT(rotationMatrix->phi);
        rotationMatrix->c = calcC(rotationMatrix->t);
        rotationMatrix->s = rotationMatrix->t * rotationMatrix->c; /* s = tc */

        /* calculate A' */
        updatedMatrix = updateAtag(matrix, rotationMatrix, rows);

        /* calculate P */
        rotMatrix = genRotationMatrix(rotationMatrix, rows);

        /* update V */
        updatedEigVectors = matrixProduct(eigenVectors, rotMatrix, rows);

        freeMatrix(eigenVectors, rows);
        eigenVectors = updatedEigVectors;

        isConverged = (iterCounter >= ROTATION_LIMIT ||
                       offDiagEuclidDistance(matrix, rows) - offDiagEuclidDistance(updatedMatrix, rows) <= EPSILON);

        if (iterCounter >= 2)
        {
            freeMatrix(matrix, rows);
        }

        freeMatrix(rotMatrix, rows);

        matrix = updatedMatrix;
    }

    eigs = matrixToEig(matrix, eigenVectors, rows);

    free(rotationMatrix);
    freeMatrix(eigenVectors, rows);
    freeMatrix(matrix, rows);

    return eigs;
}

/**
 * prints a matrix with shape = (rows, cols)
 */
void printMatrix(double **matrix, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns - 1; j++)
        {
            printf("%.4f,", matrix[i][j]);
        }
        printf("%.4f\n", matrix[i][j]);
    }
}

/**
 * prints eigenvalues and eigenvector form given eigHolder struct
 * according to the assignments instructions
 */
void printEigs(eigenHolder **eigs, int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)
    {
        printf("%.4f,", eigs[i]->eigValue);
    }

    printf("%.4f\n", eigs[i]->eigValue);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - 1; j++)
        {
            printf("%.4f,", eigs[j]->eigenVector[i]);
        }

        printf("%.4f\n", eigs[j]->eigenVector[i]);
    }
}

/**
 * simply frees array of array of doubles
 */
void freeMatrix(double **matrix, int rows)
{
    int i;

    for (i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }

    free(matrix);
}

void freeEigenHolder(eigenHolder **eigs, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        free(eigs[i]->eigenVector);
        free(eigs[i]);
    }

    free(eigs);
}