#include <Python.h>
#include "utils.h"
#include "kmeanspp.h"

double **pyListToCMatrix(PyObject *pyMatrix)
{
    int rows = PyObject_Length(pyMatrix);
    double **cMatrix = malloc(rows * sizeof(double *));

    for (int i = 0; i < rows; i++)
    {
        PyObject *row = PyList_GetItem(pyMatrix, i);
        int cols = PyObject_Length(row);

        cMatrix[i] = malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++)
        {
            cMatrix[i][j] = PyFloat_AsDouble(PyList_GetItem(row, j));
        }
    }

    return cMatrix;
}

static PyObject *cMatrixToPyList(double **pyMatrix, int rows, int columns)
{
    PyObject *pyList = PyList_New(rows);

    for (int i = 0; i < rows; i++)
    {
        PyObject *pyRow = PyList_New(columns);
        for (int j = 0; j < columns; j++)
        {
            PyList_SetItem(pyRow, j, Py_BuildValue("d", pyMatrix[i][j]));
        }

        PyList_SetItem(pyList, i, pyRow);
    }

    return pyList;
}

static PyObject *eigHolderStructToPyList(eigenHolder **eigs, int n)
{
    PyObject *pyList = PyList_New(n + 1);

    PyObject *eigenvalues = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyList_SetItem(eigenvalues, i, Py_BuildValue("d", eigs[i]->eigValue));
    }
    PyList_SetItem(pyList, 0, eigenvalues);

    for (int i = 0; i < n; i++)
    {
        PyObject *eigenvector = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            PyList_SetItem(eigenvector, j, Py_BuildValue("d", eigs[j]->eigenVector[i]));
        }

        PyList_SetItem(pyList, i + 1, eigenvector);
    }

    return pyList;
}

static PyObject *wam(PyObject *self, PyObject *args)
{
    PyObject *pyMatrix;
    double **cMatrix;

    if (!PyArg_ParseTuple(args, "O", &pyMatrix))
    {
        return NULL;
    }

    cMatrix = pyListToCMatrix(pyMatrix);
    int rows = PyObject_Length(pyMatrix);
    int cols = PyObject_Length(PyList_GetItem(pyMatrix, 0));

    return cMatrixToPyList(genWAM(cMatrix, rows, cols),
                           rows,
                           rows);
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *pyMatrix;
    double **cMatrix;

    if (!PyArg_ParseTuple(args, "O", &pyMatrix))
    {
        return NULL;
    }

    cMatrix = pyListToCMatrix(pyMatrix);
    int rows = PyObject_Length(pyMatrix);
    int cols = PyObject_Length(PyList_GetItem(pyMatrix, 0));

    return cMatrixToPyList(genDDG(
                               genWAM(cMatrix, rows, cols), rows),
                           rows,
                           rows);
}

static PyObject *gl(PyObject *self, PyObject *args)
{
    PyObject *pyMatrix;
    double **cMatrix;
    if (!PyArg_ParseTuple(args, "O", &pyMatrix))
    {
        return NULL;
    }

    cMatrix = pyListToCMatrix(pyMatrix);
    int rows = PyObject_Length(pyMatrix);
    int cols = PyObject_Length(PyList_GetItem(pyMatrix, 0));
    double **wam = genWAM(cMatrix, rows, cols);
    double **ddg = genDDG(wam, rows);

    return cMatrixToPyList(genGL(wam, ddg, rows),
                           rows,
                           rows);
}

static PyObject *jacobiHandler(PyObject *self, PyObject *args)
{
    PyObject *pyMatrix;
    double **cMatrix;

    if (!PyArg_ParseTuple(args, "O", &pyMatrix))
    {
        return NULL;
    }

    cMatrix = pyListToCMatrix(pyMatrix);
    int rows = PyObject_Length(pyMatrix);

    return eigHolderStructToPyList(jacobi(cMatrix, rows),
                                   rows);
}

static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *input_mat;
    PyObject *centroids;
    PyObject *line;
    PyObject *line_2;
    PyObject *result;
    double obj;
    double **mat = NULL;
    double **curr_clusters;
    double eps;

    int iter;
    int d, n, k, j;
    int i = 0;
    if (!PyArg_ParseTuple(args, "iidiiOO", &k, &iter, &eps, &n, &d, &input_mat, &centroids))
    {
        printf("Invalid Input! \n");
        return NULL;
    }
    if (!PyList_Check(input_mat) || !PyList_Check(centroids))
    {

        printf("Invalid Input! \n");
        return NULL;
    }
    if (k > n)
    {
        printf("Invalid Input! \n");
        return NULL;
    }
    mat = empty_matrix(n, d);
    if (mat == NULL)
    {
        printf("An error has occurred\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < n; i++)
    {
        line = PyList_GetItem(input_mat, i);
        for (j = 0; j < d; j++)
        {
            obj = PyFloat_AsDouble(PyList_GetItem(line, j));
            mat[i][j] = obj;
        }
    }
    curr_clusters = empty_matrix(k, d);
    if (curr_clusters == NULL)
    {
        printf("An error has occurred\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < k; i++)
    {
        line_2 = PyList_GetItem(centroids, i);
        for (j = 0; j < d; j++)
        {
            obj = PyFloat_AsDouble(PyList_GetItem(line_2, j));
            curr_clusters[i][j] = obj;
        }
    }
    curr_clusters = kmeans(mat, curr_clusters, n, k, d, iter, eps);
    result = PyList_New(k);
    if (result == NULL)
    {
        return NULL;
    }
    for (i = 0; i < k; i++)
    {
        line = PyList_New(d);
        if (line == NULL)
        {
            return NULL;
        }
        for (j = 0; j < d; j++)
        {
            PyList_SetItem(line, j, PyFloat_FromDouble(curr_clusters[i][j]));
        }
        PyList_SetItem(result, i, line);
    }
    freeMemory(mat, n);
    freeMemory(curr_clusters, k);
    return Py_BuildValue("O", result);
}

static PyMethodDef myMethods[] = {
    {"wam",
     (PyCFunction)wam,
     METH_VARARGS,
     PyDoc_STR("Calculates the WAM for an input matrix. input: matrix: List[List[float]] ouput: List[List[float]]")},

    {"ddg",
     (PyCFunction)ddg,
     METH_VARARGS,
     PyDoc_STR("Calculates the DDG for an input matrix. input: matrix: List[List[float]] ouput: List[List[float]]")},

    {"gl",
     (PyCFunction)gl,
     METH_VARARGS,
     PyDoc_STR("Calculates the GL for an input matrix. input: matrix: List[List[float]] ouput: List[List[float]]")},

    {"jacobi",
     (PyCFunction)jacobiHandler,
     METH_VARARGS,
     PyDoc_STR("Calculates the eigenvalues & eigenvectors for an input matrix using jacobi algorithm. input: matrix: List[List[float]] ouput: List[List[float]]")},

    {"fit",
     (PyCFunction)fit,
     METH_VARARGS,
     PyDoc_STR("A kmeans implemention written in C")},

    {NULL, NULL, 0, NULL}};

static struct PyModuleDef mykmeanssp = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    myMethods};

PyMODINIT_FUNC PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m)
    {
        return NULL;
    }
    return m;
}