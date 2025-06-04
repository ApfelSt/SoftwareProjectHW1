/*import os
import sys

def kmeans(X, k, max_iters, eps=1e-3):
    """
    Perform K-means clustering on the dataset X.

    Parameters:
    X : array-like of points in R^d
        The dataset to cluster.
    k : int
    max_iters : int
    eps : float
    """

    # The initial centroids are the first k points in X
    centroids = [X[i] for i in range(k)]
    
    for _ in range(max_iters):
        # Assign points to the nearest centroid
        clusters = [[] for _ in range(k)]
        for point in X:
            distances = [sum((p - c) ** 2 for p, c in zip(point, centroid)) for centroid in centroids]
            closest_centroid = distances.index(min(distances))
            clusters[closest_centroid].append(point)

        # Update centroids
        new_centroids = []
        for cluster in clusters:
            if cluster:  # Avoid division by zero
                new_centroid = [sum(dim) / len(cluster) for dim in zip(*cluster)] 
                new_centroids.append(new_centroid)
            else:
                new_centroids.append(random.choice(X))  # Reinitialize if empty

        # Check convergence
        if all(sum((n - o) ** 2 for n, o in zip(new_c, old_c)) < eps for new_c, old_c in zip(new_centroids, centroids)):
            centroids = new_centroids
            break
        
        centroids = new_centroids

    return centroids


def main():
    """Usage: python3 kmeans.py <k> <max_iters> < <data_file>"""
    if len(sys.argv) != 3:
        print(main.__doc__)
        sys.exit(1)

    k = int(sys.argv[1])
    max_iters = int(sys.argv[2])

    # Read data from stdin
    X = []
    for line in sys.stdin:
        # split by a comma
        point = list(map(float, line.strip().split(',')))
        X.append(point)

    # Perform K-means clustering
    centroids = kmeans(X, k, max_iters)

    # Print the final centroids with 4 decimal places
    for centroid in centroids:
        print(','.join(f"{c:.4f}" for c in centroid))


if __name__ == "__main__":
    main()*/

/* that was previously provided is a Python implementation of the K-means clustering algorithm. Below is the equivalent C code that performs the same functionality:
 */

/*---------------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Code
void kmeans(double **X, int n, int d, int k, int max_iters, double eps, double **centroids) {
    // Initialize centroids with the first k points
    for (int i = 0; i < k; i++) {
        memcpy(centroids[i], X[i], d * sizeof(double));
    }

    for (int iter = 0; iter < max_iters; iter++) {
        // Assign points to the nearest centroid
        int *clusters = malloc(n * sizeof(int));
        for (int i = 0; i < n; i++) {
            double min_dist = INFINITY;
            for (int j = 0; j < k; j++) {
                double dist = 0.0;
                for (int l = 0; l < d; l++) {
                    dist += (X[i][l] - centroids[j][l]) * (X[i][l] - centroids[j][l]);
                }
                if (dist < min_dist) {
                    min_dist = dist;
                    clusters[i] = j;
                }
            }
        }

        // Update centroids
        double **new_centroids = malloc(k * sizeof(double *));
        for (int j = 0; j < k; j++) {
            new_centroids[j] = calloc(d, sizeof(double));
        }
        int *counts = calloc(k, sizeof(int));

        for (int i = 0; i < n; i++) {
            int cluster_id = clusters[i];
            counts[cluster_id]++;
            for (int l = 0; l < d; l++) {
                new_centroids[cluster_id][l] += X[i][l];
            }
        }

        // Average the new centroids
        for (int j = 0; j < k; j++) {
            if (counts[j] > 0) {
                for (int l = 0; l < d; l++) {
                    new_centroids[j][l] /= counts[j];
                }
            } else {
                // Reinitialize if empty
                for (int l = 0; l < d; l++) {
                    new_centroids[j][l] = X[rand() % n][l];
                }
            }
        }

        // Check convergence
        int converged = 1;
        for (int j = 0; j < k; j++) {
            double dist = 0.0;
            for (int l = 0; l < d; l++) {
                dist += (new_centroids[j][l] - centroids[j][l]) * (new_centroids[j][l] - centroids[j][l]);
            }
            if (dist >= eps * eps) {
                converged = 0;
            }
        }
        if (converged) {
            for (int j = 0; j < k; j++) {
                memcpy(centroids[j], new_centroids[j], d * sizeof(double));
            }
            for (int j = 0; j < k; j++) {
                free(new_centroids[j]);
            }
            free(new_centroids);
            free(counts);
            free(clusters);
            return;
        }

        // Update centroids
        for (int j = 0; j < k; j++) {
            memcpy(centroids[j], new_centroids[j], d * sizeof(double));
        }
        for (int j = 0; j < k; j++) {
            free(new_centroids[j]);
        }
        free(new_centroids);
        free(counts);
        free(clusters);
    }
    return;
}

void print_centroids(double **centroids, int k, int d) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < d; j++) {
            printf("%.4f", centroids[i][j]);
            if (j < d - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <k> <max_iters>\n", argv[0]);
        return 1;
    }

    int k = atoi(argv[1]);
    int max_iters = atoi(argv[2]);

    // Read data from stdin
    double **X = malloc(1000 * sizeof(double *)); // Assuming a maximum of 1000 points
    int n = 0, d = 0;

    char line[1024];
    while (fgets(line, sizeof(line), stdin)) {
        X[n] = malloc(10 * sizeof(double)); // Assuming a maximum of 10 dimensions
        char *token = strtok(line, ",");
        int dim = 0;
        while (token) {
            X[n][dim++] = atof(token);
            token = strtok(NULL, ",");
        }
        if (d == 0) {
            d = dim; // Set dimension on first read
        }
        n++;
    }

    // Allocate memory for centroids
    double **centroids = malloc(k * sizeof(double *));
    for (int i = 0; i < k; i++) {
        centroids[i] = malloc(d * sizeof(double));
    }

    // Perform K-means clustering
    kmeans(X, n, d, k, max_iters, 1e-3, centroids);
    print_centroids(centroids, k, d);

    // Free allocated memory
    for (int i = 0; i < n; i++) {
        free(X[i]);
    }
    free(X);
    
    for (int i = 0; i < k; i++) {
        free(centroids[i]);
    }
    free(centroids);

    return 0;
}

----------------------------------------------------------------
*/
/* TODO: check that we do not use float and if so then use double instead */
/* TODO: only use stdlib.h, math.h, and stdio.h. currently we use string.h for memcpy */
/* TODO: check any other requirements for the C code */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#ifndef INFINITY
#define INFINITY DBL_MAX
#endif

void kmeans(double **X, int n, int d, int k, int max_iters, double eps, double **centroids) {
    int *clusters;
    double **new_centroids;
    int *counts;
    int i, j, l, iter, cluster_id, converged;
    double min_dist, dist;

    clusters = malloc(n * sizeof(int));
    if (clusters == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (i = 0; i < k; i++) {
        memcpy(centroids[i], X[i], d * sizeof(double));
    }

    for (iter = 0; iter < max_iters; iter++) {
        for (i = 0; i < n; i++) {
            min_dist = INFINITY;
            for (j = 0; j < k; j++) {
                dist = 0.0;
                for (l = 0; l < d; l++) {
                    dist += (X[i][l] - centroids[j][l]) * (X[i][l] - centroids[j][l]);
                }
                if (dist < min_dist) {
                    min_dist = dist;
                    clusters[i] = j;
                }
            }
        }

        new_centroids = malloc(k * sizeof(double *));
        if (new_centroids == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            free(clusters);
            exit(1);
        }
        for (j = 0; j < k; j++) {
            new_centroids[j] = calloc(d, sizeof(double));
            if (new_centroids[j] == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                for (i = 0; i < j; i++) {
                    free(new_centroids[i]);
                }
                free(new_centroids);
                free(clusters);
                exit(1);
            }
        }

        counts = calloc(k, sizeof(int));
        if (counts == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            for (j = 0; j < k; j++) {
                free(new_centroids[j]);
            }
            free(new_centroids);
            free(clusters);
            exit(1);
        }

        for (i = 0; i < n; i++) {
            cluster_id = clusters[i];
            counts[cluster_id]++;
            for (l = 0; l < d; l++) {
                new_centroids[cluster_id][l] += X[i][l];
            }
        }

        for (j = 0; j < k; j++) {
            if (counts[j] > 0) {
                for (l = 0; l < d; l++) {
                    new_centroids[j][l] /= counts[j];
                }
            } else {
                for (l = 0; l < d; l++) {
                    new_centroids[j][l] = X[rand() % n][l];
                }
            }
        }

        converged = 1;
        for (j = 0; j < k; j++) {
            dist = 0.0;
            for (l = 0; l < d; l++) {
                dist += (new_centroids[j][l] - centroids[j][l]) * (new_centroids[j][l] - centroids[j][l]);
            }
            if (dist >= eps * eps) {
                converged = 0;
            }
        }

        for (j = 0; j < k; j++) {
            memcpy(centroids[j], new_centroids[j], d * sizeof(double));
        }

        for (j = 0; j < k; j++) {
            free(new_centroids[j]);
        }
        free(new_centroids);
        free(counts);

        if (converged) {
            break;
        }
    }

    free(clusters);
}

void print_centroids(double **centroids, int k, int d) {
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            printf("%.4f", centroids[i][j]);
            if (j < d - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

/* int main(int argc, char *argv[]) {
    int k, max_iters, n, d, i, dim;
    double **X, **centroids;
    char line[1024];
    char *token;

    if (argc != 3 && argc != 2) {
        fprintf(stderr, "Usage: %s <k> <max_iters> or %s <k>\n", argv[0], argv[0]);
        return 1;
    }

    k = atoi(argv[1]);
    max_iters = (argc == 3) ? atoi(argv[2]) : 400;

    n = 0;
    d = 0;
    while (fgets(line, sizeof(line), stdin)) {
        n++;
        if (n == 1) {
            token = strtok(line, ",");
            while (token) {
                d++;
                token = strtok(NULL, ",");
            }
        }
    }

    X = malloc(n * sizeof(double *));
    if (X == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    n = 0;
    d = 0;
    while (fgets(line, sizeof(line), stdin)) {
        X[n] = malloc(d * sizeof(double));
        if (X[n] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            for (i = 0; i < n; i++) {
                free(X[i]);
            }
            free(X);
            return 1;
        }
        token = strtok(line, ",");
        dim = 0;
        while (token) {
            X[n][dim++] = atof(token);
            token = strtok(NULL, ",");
        }
        if (d == 0) {
            d = dim;
        }
        n++;
    }

    centroids = malloc(k * sizeof(double *));
    if (centroids == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        for (i = 0; i < n; i++) {
            free(X[i]);
        }
        free(X);
        return 1;
    }
    for (i = 0; i < k; i++) {
        centroids[i] = malloc(d * sizeof(double));
        if (centroids[i] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            for (dim = 0; dim < i; dim++) {
                free(centroids[dim]);
            }
            free(centroids);
            for (dim = 0; dim < n; dim++) {
                free(X[dim]);
            }
            free(X);
            return 1;
        }
    }

    kmeans(X, n, d, k, max_iters, 1e-3, centroids);
    print_centroids(centroids, k, d);

    for (i = 0; i < n; i++) {
        free(X[i]);
    }
    free(X);

    for (i = 0; i < k; i++) {
        free(centroids[i]);
    }
    free(centroids);

    return 0;
} */

int main(int argc, char *argv[]) {
    int k, max_iters, n, d, i, dim, capacity;
    double **X, **centroids;
    char line[1024];
    char *token;

    if (argc != 3 && argc != 2) {
        fprintf(stderr, "Usage: %s <k> <max_iters> or %s <k>\n", argv[0], argv[0]);
        return 1;
    }

    k = atoi(argv[1]);
    max_iters = (argc == 3) ? atoi(argv[2]) : 400;

    /* TODO: Check if k and max_iters are within valid ranges */

    /* Initialize variables */
    n = 0;
    d = 0;
    capacity = 10;
    X = malloc(capacity * sizeof(double *));
    if (X == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    /* Read points from stdin */
    while (fgets(line, sizeof(line), stdin)) {
        if (n >= capacity) {
            capacity *= 2;
            X = realloc(X, capacity * sizeof(double *));
            if (X == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                for (i = 0; i < n; i++) {
                    free(X[i]);
                }
                free(X);
                return 1;
            }
        }

        X[n] = NULL;
        token = strtok(line, ",");
        dim = 0;

        /* Parse the line and allocate memory for dimensions */
        while (token) {
            if (dim == 0) {
                X[n] = malloc(10 * sizeof(double));
                if (X[n] == NULL) {
                    fprintf(stderr, "Memory allocation failed\n");
                    for (i = 0; i < n; i++) {
                        free(X[i]);
                    }
                    free(X);
                    return 1;
                }
            } else if (dim % 10 == 0) {
                X[n] = realloc(X[n], (dim + 10) * sizeof(double));
                if (X[n] == NULL) {
                    fprintf(stderr, "Memory allocation failed\n");
                    for (i = 0; i < n; i++) {
                        free(X[i]);
                    }
                    free(X);
                    return 1;
                }
            }
            X[n][dim++] = atof(token);
            token = strtok(NULL, ",");
        }

        if (d == 0) {
            d = dim; 
        } else if (dim != d) {
            fprintf(stderr, "Inconsistent dimensions in input\n");
            for (i = 0; i <= n; i++) {
                free(X[i]);
            }
            free(X);
            return 1;
        }

        n++;
    }

    /* Allocate memory for centroids */
    centroids = malloc(k * sizeof(double *));
    if (centroids == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        for (i = 0; i < n; i++) {
            free(X[i]);
        }
        free(X);
        return 1;
    }
    for (i = 0; i < k; i++) {
        centroids[i] = malloc(d * sizeof(double));
        if (centroids[i] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            for (dim = 0; dim < i; dim++) {
                free(centroids[dim]);
            }
            free(centroids);
            for (dim = 0; dim < n; dim++) {
                free(X[dim]);
            }
            free(X);
            return 1;
        }
    }

    /* Perform K-means clustering */
    kmeans(X, n, d, k, max_iters, 1e-3, centroids);
    print_centroids(centroids, k, d);

    /* Free allocated memory */
    for (i = 0; i < n; i++) {
        free(X[i]);
    }
    free(X);

    for (i = 0; i < k; i++) {
        free(centroids[i]);
    }
    free(centroids);

    return 0;
}
