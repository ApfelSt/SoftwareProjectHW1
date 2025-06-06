#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int integer_check_in_range (const char *str, int lower_limit, int higer_limit) {
    int value;
    value = atoi(str);
    if (value <= lower_limit || value >= higer_limit) {
        return 0;
    }
    if (*str == '\0') {
        return 0;
    }
    while (*str && *str != '.') {
        if (*str < '0' || *str > '9') {
            return 0;
        }
        str++;
    }
    if (*str == '.') {
        while (*str) {
            if (*str != '0') {
                return 0;
            }
            str++;
        }
    }
    return 1;
}

void kmeans(double **X, int n, int d, int k, int max_iters, double eps, double **centroids) {
    int *clusters;
    double **new_centroids;
    int *counts;
    int i, j, l, iter, cluster_id, converged;
    double min_dist, dist;

    clusters = malloc(n * sizeof(int));
    if (clusters == NULL) {
        fprintf(stderr, "An Error Has Occurred\n");
        exit(1);
    }

    for (i = 0; i < k; i++) {
        memcpy(centroids[i], X[i], d * sizeof(double));
    }

    for (iter = 0; iter < max_iters; iter++) {
        for (i = 0; i < n; i++) {
            min_dist = -1.0;
            for (j = 0; j < k; j++) {
                dist = 0.0;
                for (l = 0; l < d; l++) {
                    dist += (X[i][l] - centroids[j][l]) * (X[i][l] - centroids[j][l]);
                }
                if (dist < min_dist || min_dist < 0.0) {
                    min_dist = dist;
                    clusters[i] = j;
                }
            }
        }

        new_centroids = malloc(k * sizeof(double *));
        if (new_centroids == NULL) {
            fprintf(stderr, "An Error Has Occurred\n");
            free(clusters);
            exit(1);
        }
        for (j = 0; j < k; j++) {
            new_centroids[j] = calloc(d, sizeof(double));
            if (new_centroids[j] == NULL) {
                fprintf(stderr, "An Error Has Occurred\n");
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
            fprintf(stderr, "An Error Has Occurred\n");
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


    /* Initialize variables */
    n = 0;
    d = 0;
    capacity = 10;
    X = malloc(capacity * sizeof(double *));
    if (X == NULL) {
        fprintf(stderr, "An Error Has Occurred\n");
        return 1;
    }

    /* Read points from stdin */
    while (fgets(line, sizeof(line), stdin)) {
        if (n >= capacity) {
            capacity *= 2;
            X = realloc(X, capacity * sizeof(double *));
            if (X == NULL) {
                fprintf(stderr, "An Error Has Occurred\n");
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
                    fprintf(stderr, "An Error Has Occurred\n");
                    for (i = 0; i < n; i++) {
                        free(X[i]);
                    }
                    free(X);
                    return 1;
                }
            } else if (dim % 10 == 0) {
                X[n] = realloc(X[n], (dim + 10) * sizeof(double));
                if (X[n] == NULL) {
                    fprintf(stderr, "An Error Has Occurred\n");
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
            fprintf(stderr, "An Error Has Occurred\n");
            for (i = 0; i <= n; i++) {
                free(X[i]);
            }
            free(X);
            return 1;
        }

        n++;
    }

    /*Checking input*/
    if (!integer_check_in_range(argv[1], 1, n)) {
        printf("Incorrect number of clusters!\n");
        for (i = 0; i < n; i++) {
            free(X[i]);
        }
        free(X);
        return 1;
    }
    
    if (argc == 3 && !integer_check_in_range(argv[2], 1, 1000)) {
        printf("Incorrect maximum iteration!\n");
        for (i = 0; i < n; i++) {
            free(X[i]);
        }
        free(X);
        return 1;
    }

    /* Allocate memory for centroids */
    centroids = malloc(k * sizeof(double *));
    if (centroids == NULL) {
        fprintf(stderr, "An Error Has Occurred\n");
        for (i = 0; i < n; i++) {
            free(X[i]);
        }
        free(X);
        return 1;
    }
    for (i = 0; i < k; i++) {
        centroids[i] = malloc(d * sizeof(double));
        if (centroids[i] == NULL) {
            fprintf(stderr, "An Error Has Occurred\n");
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
