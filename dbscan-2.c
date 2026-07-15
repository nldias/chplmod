#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define UNVISITED 0
#define NOISE -1

typedef struct {
    double x, y;
} Point;

double euclidean_distance(Point p1, Point p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Find all points within epsilon of a given point
void find_neighbors(int p_idx, const Point* points, int num_points, double eps, int* neighbors, int* num_neighbors) {
    *num_neighbors = 0;
    for (int i = 0; i < num_points; ++i) {
        if (p_idx == i) continue;
        if (euclidean_distance(points[p_idx], points[i]) <= eps) {
            neighbors[(*num_neighbors)++] = i;
        }
    }
}

// Expand a cluster from a core point
void expand_cluster(
   int p_idx,
   int cluster_id,
   Point* points,
   int num_points,
   double eps,
   int minPts,
   int* labels,
   int* neighbors_buffer,
   int* queue_buffer) {
    int queue_head = 0, queue_tail = 0;
    queue_buffer[queue_tail++] = p_idx;
    labels[p_idx] = cluster_id;

    while (queue_head < queue_tail) {
        int current_idx = queue_buffer[queue_head++];
        int num_current_neighbors = 0;
        find_neighbors(current_idx, points, num_points, eps, neighbors_buffer, &num_current_neighbors);

        if (num_current_neighbors >= minPts - 1) { // -1 because find_neighbors excludes itself
            for (int i = 0; i < num_current_neighbors; ++i) {
                int neighbor_idx = neighbors_buffer[i];
                if (labels[neighbor_idx] == UNVISITED) {
                    labels[neighbor_idx] = cluster_id;
                    queue_buffer[queue_tail++] = neighbor_idx;
                } else if (labels[neighbor_idx] == NOISE) {
                    labels[neighbor_idx] = cluster_id; // Change noise to border point
                }
            }
        }
    }
}

void dbscan(Point* points, int num_points, double eps, int minPts, int* labels) {
    int cluster_id = 0;
    // Buffers for neighbors and queue
    int* neighbors_buffer = (int*)malloc(num_points * sizeof(int));
    int* queue_buffer = (int*)malloc(num_points * sizeof(int));

    for (int i = 0; i < num_points; ++i) {
        if (labels[i] != UNVISITED) {
            continue; // Already processed
        }

        int num_neighbors = 0;
        find_neighbors(i, points, num_points, eps, neighbors_buffer, &num_neighbors);

        if (num_neighbors < minPts - 1) { // If fewer than minPts neighbors (excluding self)
            labels[i] = NOISE;
        } else {
            cluster_id++;
            expand_cluster(i, cluster_id, points, num_points, eps, minPts, labels, neighbors_buffer, queue_buffer);
        }
    }

    free(neighbors_buffer);
    free(queue_buffer);
}

// Example Usage
int main() {
    Point data[] = {
        {0.1, 1.0}, {0.2, 0.9}, {0.3, 1.0}, {0.4, 0.6}, {0.5, 0.6},
        {0.6, 0.5}, {0.7, 0.8}, {0.8, 0.1}, {0.9, 0.2}, {1.0, 0.1}
    };
    int num_points = sizeof(data) / sizeof(data[0]);
    int labels[num_points];
    // Initialize labels to UNVISITED (0)
    for (int i = 0; i < num_points; ++i) labels[i] = UNVISITED;

    double epsilon = 0.20;
    int minPts = 2; // DBSCAN parameters

    printf("Clustering with epsilon = %.2f and minPts = %d\n", epsilon, minPts);
    dbscan(data, num_points, epsilon, minPts, labels);

    printf("\nClustering results:\n");
    for (int i = 0; i < num_points; ++i) {
        printf("Point (%.1f, %.1f) | Cluster ID: %d\n", data[i].x, data[i].y, labels[i]);
    }

    return 0;
}
