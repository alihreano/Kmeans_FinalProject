import spkmeansmodule as km
import sys
import pandas as pd
import numpy as np


def kmeans_plusplus(data, k):
    """
    Perform k-means++ initialization
    :param data: numpy array of data points
    :param k: number of clusters
    :return: numpy array of initialized centroids and inidces of the points chosen for the inital clusters
    """
    np.random.seed(0)
    n, d = data.shape

    index_arr = np.empty((k))
    first_centroid = np.random.choice(n)
    index_arr[0] = int(first_centroid)
    centroids = np.empty((k, d))
    centroids[0] = data[first_centroid]
    distances = np.empty(n)

    for i in range(1, k):
        # for each data point, find the distance to the closest centroid
        for j, point in enumerate(data):
            tmp = point[:]
            # tmp = point[1:]
            distances[j] = np.min([np.linalg.norm(tmp - c)
                                  for c in centroids[:i]])
        # normalize the distance array
        distances /= np.sum(distances)
        # choose next centroid with probability proportional to distance
        cumulative = np.cumsum(distances)
        rand_num = np.random.rand()
        for j, p in enumerate(cumulative):
            if rand_num <= p:
                # centroids[i] = data[j][1:]
                centroids[i] = data[j][:]
                index_arr[i] = j
                break

    return centroids, index_arr


def execute(k, maxItr, epsilon, input_matrix):
    """calculate centroids with kmeansand pribnts the indices and final centroids

    Args:
        k (int): number of clusters.
        maxItr (int): max iteration we want the kmeans to run
        epsilon (float): epsilon for convergence in kmeans
        input_matrix (np.array): matrix of first k eigenvalues calculated from the original input

    returns: void
    """
    input_array = input_matrix.tolist()
    mat = []

    for point in input_array:
        mat.append(point[:])

    n = len(mat)
    d = len(mat[0])

    if (k >= n) or (k <= 1) or (int(k) != k):
        print("Invalid number of clusters!")
        return
    if (maxItr <= 1) or (maxItr >= 1000) or (int(maxItr) != maxItr):
        print("Invalid maximum iteration!")
        return

    # centroids µ1, µ2, ... , µK ∈ R^d where 1<K<N.
    centroids, index_arr = kmeans_plusplus(input_matrix, k)
    matrix = km.fit(
        int(k), maxItr, epsilon, n, d, mat, centroids.tolist())
    return index_arr, matrix


def eigengap_heuristic(eigen_values):
    mid = len(eigen_values) // 2  # int divison
    gaps = np.zeros(mid)
    gaps = np.array([np.abs(eigen_values[i] - eigen_values[i + 1])
                    for i in range(mid)])
    return np.argmax(gaps) + 1


def spkmeans(matrix, k):
    """pre-processing for kmeans++: calc and sort eigenvalue, takes only
    first k. (if k == -1 we use eigengap heuristic)

    Args:
        matrix (List[List[float]]): input data
        k (int): num of clusters. can be -1 if k not provided
    """
    gl = km.gl(matrix)
    eigs = np.array(km.jacobi(gl))

    # sort the eigen vectors according to the eigenvalues
    sorted_eigs = eigs[:, eigs[0].argsort()]

    if k == -1:
        k = eigengap_heuristic(sorted_eigs[0])

    filtered_eigenvectors = pd.DataFrame(sorted_eigs[1:, :k]).values
    return execute(k, 300, 0.0, filtered_eigenvectors)


def main(argc, argv):
    def print_vector(vector): return print(
        ",".join([str(int(entry)) for entry in vector]))

    def print_matrix(matrix): return print(
        "\n".join([",".join([("%.4f" % entry) for entry in row]) for row in matrix]))

    if argc == 4:
        k, goal, file_name = argv[1:]
    elif argc == 3:
        k = -1
        goal, file_name = argv[1:]
    else:
        print("An Error Has Occurred")

    matrix = pd.read_csv(file_name, header=None).to_numpy().tolist()

    if goal == "spk":
        index_arr, result_matrix = spkmeans(matrix, k)
        print_vector(index_arr)
    elif goal == "wam":
        result_matrix = km.wam(matrix)
    elif goal == "ddg":
        result_matrix = km.ddg(matrix)
    elif goal == "gl":
        result_matrix = km.gl(matrix)
    elif goal == "jacobi":
        result_matrix = km.jacobi(matrix)

    print_matrix(result_matrix)


if __name__ == '__main__':
    main(len(sys.argv), sys.argv)
