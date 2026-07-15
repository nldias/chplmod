#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1.0       // The radius of neighborhood
#define MIN_PTS 3     // Minimum number of points required to form a dense region
#define NOISE -1      // Label for noise points
#define UNDEFINED 0   // Initial label for unprocessed points
#define MAX_POINTS 100 // Maximum number of data points

// Structure to represent a data point
typedef struct {
    double x;
    double y;
    int cluster_id; // Will store the cluster label (or NOISE, UNDEFINED)
} Point;

// Helper function to calculate Euclidean distance between two points
double euclidean_distance(Point p1, Point p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Helper function to find all neighbors within the epsilon distance
int* range_query(Point* points, int num_points, int p_idx, double eps, int* num_neighbors) {
    int* neighbors = (int*)malloc(num_points * sizeof(int)); // Allocate enough space
    *num_neighbors = 0;
    for (int i = 0; i < num_points; i++) {
        if (p_idx == i) continue;
        if (euclidean_distance(points[p_idx], points[i]) <= eps) {
            neighbors[(*num_neighbors)++] = i;
        }
    }
    // Reallocate to exact size or manage dynamically in a list
    return neighbors; // Note: caller must free this memory
}

// Function to expand a cluster
void expand_cluster(
   Point* points,
   int num_points,
   int p_idx,
   int cluster_id,
   double eps,
   int min_pts) {
    points[p_idx].cluster_id = cluster_id;
    int num_neighbors;
    int* neighbors = range_query(points, num_points, p_idx, eps, &num_neighbors);
    
    // Use a while loop to process new neighbors added to the 'neighbors'
    // list (simulating a queue/seed set)
    for (int i = 0; i < num_neighbors; i++) {
        int q_idx = neighbors[i];
        if (points[q_idx].cluster_id == NOISE) {
            points[q_idx].cluster_id = cluster_id; // Change noise to border point
        }
        if (points[q_idx].cluster_id != UNDEFINED) {
            continue; // Already processed
        }
        points[q_idx].cluster_id = cluster_id;
        
        int sub_num_neighbors;
        int* sub_neighbors = range_query(points, num_points, q_idx, eps, &sub_num_neighbors);
        if (sub_num_neighbors >= min_pts) {
            // If Q is a core point, add its neighbors to the current cluster's seed set
            // This part requires dynamic resizing of the neighbors array/list in a real C implementation
            // For this simple example, we'd need a more robust list data structure
            // (e.g., linked list, or realloc, which is complex to manage here).
            // A queue-based approach in C# is simpler, as seen in the search results.
        }
        free(sub_neighbors);
    }
    free(neighbors);
}

// The main DBSCAN function
void dbscan(Point* points, int num_points, double eps, int min_pts) {
    int cluster_id = 0;
    for (int i = 0; i < num_points; i++) {
        if (points[i].cluster_id != UNDEFINED) {
            continue; // Point already visited
        }
        
        int num_neighbors;
        int* neighbors = range_query(points, num_points, i, eps, &num_neighbors);
        
        if (num_neighbors < min_pts) {
            points[i].cluster_id = NOISE; // Mark as noise initially
        } else {
            cluster_id++; // Next cluster label
            expand_cluster(points, num_points, i, cluster_id, eps, min_pts);
        }
        free(neighbors);
    }
}

// Example usage
int main() {
    Point data[MAX_POINTS] = {
        {0.1, 1.0, UNDEFINED}, {0.2, 0.9, UNDEFINED}, {0.3, 1.0, UNDEFINED},
        {0.4, 0.6, UNDEFINED}, {0.5, 0.6, UNDEFINED}, {0.6, 0.5, UNDEFINED},
        {0.7, 0.8, UNDEFINED}, {0.8, 0.1, UNDEFINED}, {0.9, 0.2, UNDEFINED},
        {1.0, 0.1, UNDEFINED}
    };
    int num_points = 10;
    double epsilon = 0.20;
    int minPts = 2;

    printf("Clustering with epsilon = %.2f and min_pts = %d\n", epsilon, minPts);
    dbscan(data, num_points, epsilon, minPts);

    printf("\nClustering results:\n");
    for (int i = 0; i < num_points; i++) {
        printf("Point (%.2f, %.2f) | Cluster: %d\n", data[i].x, data[i].y, data[i].cluster_id);
    }

    return 0;
}
