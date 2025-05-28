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
            break
        
        centroids = new_centroids

    return centroids
