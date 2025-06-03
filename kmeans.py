import os
import sys

def kmeans(X, k, max_iters=400, eps=1e-3):
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
    """Usage: python3 kmeans.py <k> <max_iters> < <data_file> or python3 kmeans.py <k> < data_file.txt"""
    if len(sys.argv) != 3 and len(sys.argv) != 2:
        print(main.__doc__)
        sys.exit(1)

    k = int(sys.argv[1])
    if len(sys.argv) == 2:
        max_iters = 400  # Default value
    else:
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
    main()

