private const EPS = 1.0;           // the radius of neighborhood
private const MIN_PTS = 3;         // minimum number of points required to form a
                                   // dense region
private const NOISE = -1;          // label for noise points
private const UNDEFINED = 0;       // initial label for unprocessed points
private const MAX_POINTS = 100;    // maximum number of data points


record Point {                // structure to represent a data point
   var x: real;
   var y: real;
   var cluster_id: int;       // will store the cluster label (or NOISE, UNDEFINED)
} 

// -----------------------------------------------------------------------------
// --> euclidean_distance: Euclidean distance between two points
// -----------------------------------------------------------------------------
private inline proc euclidean_distance(
   const ref p1: Point,
   const ref p2: Point
   ): real {
    return sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2);
}

// -----------------------------------------------------------------------------
// --> range_query: Helper function to find all neighbors within the epsilon
// distance of *this* point.
// -----------------------------------------------------------------------------
private proc range_query(
   const in p_idx: int,            // index of *this* point
   const in eps: real,             // radius of search
   const ref points: [?dpts] Point // array with all points
   ): [] int where points.rank == 1 {
   var dnei = {1..10};             // domain of neighborhood
   var neighbors = [dnei] int;     // neighborhood (Allocate enough space)
   var num_neighbors = 0;          // # of neigbors
   for i in dpts do {              // Loop over all points.
      if (p_idx == i) then continue;    // Skip *this* point.
      if (euclidean_distance(points[p_idx], points[i]) <= eps) then {
         num_neighbors += 1;            // Increment # of neighbors.
         dgrow(num_neighbors,dnei);     // Adjust `dnei` if necessary.
         neighbors[num_neighbors] = i;  // Include `points[i]` in neighborhodd.
      }
   }
   dnei = {1..num_neighbors};
   return neighbors; // Note: caller must free this memory
}



// -----------------------------------------------------------------------------
// --> expand_cluster: Function to expand a cluster.
// -----------------------------------------------------------------------------
private proc expand_cluster(
   const in p_idx: int,                 // *this* point
   const in cluster_id: int,            // *this* cluster
   const in min_pts: int,               // min # of points in cluster
   const in eps: real,                  // radius of neighborhood
   const ref points: [?dpts] Point      // array with all points
   ) where points.rank == 1 {
   // --------------------------------------------------------------------------
   // Insert *this* point in *this* cluster.
   // --------------------------------------------------------------------------
   points[p_idx].cluster_id = cluster_id;
   const in p_idx: int,            // index of *this* point
   const in eps: real,             // radius of search
   const ref points: [?dpts] Point // array with all points
   var neighbors = range_query(p_idx, eps, points);
   ref dnei = neighbors.indices;   // neighbors's domain
   // --------------------------------------------------------------------------
   // Use a while loop to process new neighbors added to the 'neighbors' list
   // (simulating a queue/seed set).
   // --------------------------------------------------------------------------
   for i in dnei do {
      var q_idx = neighbors[i];
      if (points[q_idx].cluster_id == NOISE) then {
         points[q_idx].cluster_id = cluster_id;   // Change noise to border pt.
      }
      if (points[q_idx].cluster_id != UNDEFINED) {
         continue;                                // Already processed
      }
      // -----------------------------------------------------------------------
      // Ith point (naturally?) belongs to neighborhood.
      // -----------------------------------------------------------------------
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


/*

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
*/