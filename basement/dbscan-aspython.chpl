// 
// This is a simple implementation of DBSCAN intended to explain the algorithm.
// 
//
// @author: Chris McCormick
// 
record Point {
   var x: real;
   var y: real;
}
// -----------------------------------------------------------------------------
// --> dbscan: Cluster the dataset `D` using the DBSCAN algorithm.
// 
// dbscam takes a dataset `D` (a list of vectors), a threshold distance `eps`,
// and a required number of points `MinPts`.
//
// It will return an array of cluster labels. The label -1 means noise, and then
// the clusters are numbered starting from 1.
// -----------------------------------------------------------------------------
proc dbscan(
   const in MinPts: int,      // The minimum # of pts in a cluster
   const in eps: real,        // The radius of a cluster
   const ref D: [?dpts] Point // The data points and their domain dpts
   ): [dpts] int where D.rank == 1 {
   // -------------------------------------------------------------------------- 
   // This list will hold the final cluster assignment for each point in D.
   // There are two reserved values:
   //    -1 - Indicates a noise point
   //     0 - Means the point hasn't been considered yet.
   // Initially all labels are 0.    
   // --------------------------------------------------------------------------
   var ptlabel: [dpts] int = 0;    // the cluster id of each point (initially 0)
   var C = 0;                      // C is the ID of the current cluster.
   // --------------------------------------------------------------------------
   // This outer loop is just responsible for picking new seed points--a point
   // from which to grow a new cluster.  Once a valid seed point is found, a new
   // cluster is created, and the cluster growth is all handled by the
   // `expandCluster` routine.
   //    
   // For each point with index iP in the Dataset D, do: 
   // (`iP` is the index of the datapoint, rather than the datapoint itself.)
   // --------------------------------------------------------------------------
   for iP in dpts do {
      // -----------------------------------------------------------------------
      // Only points that have not already been claimed can be picked as new
      // seed points.  If the point's label is not 0, continue to the next
      // point.
      // -----------------------------------------------------------------------
      if !(ptlabels[iP] == 0) then {
         continue;
      }
      // -----------------------------------------------------------------------
      // Find all of P's neighboring points.
      // -----------------------------------------------------------------------
      var NeighborPts = regionQuery(iP, eps, D);
      // -----------------------------------------------------------------------
      // If the number is below MinPts, this point is noise.  This is the only
      // condition under which a point is labeled NOISE --- when it's not a
      // valid seed point. A NOISE point may later be picked up by another
      // cluster as a boundary point (this is the only condition under which a
      // cluster label can change --- from NOISE to something else).
      // -----------------------------------------------------------------------
      if NeighborPts.size < MinPts then {
         ptlabel[iP] = -1;
      }
      // -----------------------------------------------------------------------
      // Otherwise, if there are at least MinPts nearby, use this point as the 
      // seed for a new cluster.    
      // -----------------------------------------------------------------------
      else {
         C += 1;
         growCluster(iP, C, MinPts, eps, D, ptlabel, NeighborPts);
      }
      // -----------------------------------------------------------------------
      // All data has been clustered!
      // -----------------------------------------------------------------------
   }
   return ptlabel;
}
// -----------------------------------------------------------------------------
// --> regionQuery: Find all points in dataset `D` within distance `eps` of
// point `P`.
//
// This function calculates the distance between a point P and every other 
// point in the dataset, and then returns only those points which are within a
// threshold distance `eps`.
// -----------------------------------------------------------------------------
proc regionQuery(
   const in iP: int,
   const in eps: real,
   const ref D: [dpts] Point
   ): [] int where ( D.rank == 1 ) {
   var dnei = {0..9};
   var neighbors: [dnei] int = {} ;
   // --------------------------------------------------------------------------
   //  For each point in the dataset...
   // --------------------------------------------------------------------------
   var k = -1;
   for nP in dpts do {
      // -----------------------------------------------------------------------        
      // If the distance is below the threshold, add it to the neighbors list.
      // -----------------------------------------------------------------------
      if dist(D[iP],D[nP]) < eps then {
         k += 1;
         dgrow(k,dnei);
         neighbors[k] = nP;
      }
   }
   dnei = {1..k};
   return neighbors;
}


// -----------------------------------------------------------------------------
// --> growCluster: Grow a new cluster with label `C` from the seed point `iP`.
//    
// This function searches through the dataset to find all points that belong
// to this new cluster. When this function returns, cluster `C` is complete.
// -----------------------------------------------------------------------------
proc growCluster(
   const in iP: int,               // indx of the seed pt for this new cluster
   const in C: int,                // the label for this new cluster.  
   const in MinPts: int,           // minimum required number of neighbors
   const in eps: real,             // threshold distance
   const ref D: [?dpts] Point,     // the dataset (a list of `Point`s)
   ref ptlabel: [dpts] int,        // list storing the cluster labels for all
                                   // dataset points
   ref NeighborPts [?dnei] int     // all of the neighbors of `iP`
   ) where (D.rank == 1 && ptlabel.rank == 1 && NeighborPts.rank == 1) {
   // --------------------------------------------------------------------------
   // Assert that D and NeighborPts are 0-based arrays.
   // --------------------------------------------------------------------------
   assert (dpts.low == 0 && dnei.low == 0);
   // --------------------------------------------------------------------------
   //  Assign the cluster label to the seed point.
   // --------------------------------------------------------------------------
   ptlabel[iP] = C;
   // --------------------------------------------------------------------------
   // Look at each neighbor of P (neighbors are referred to as Pn).  NeighborPts
   // will be used as set of points to search --- that is, it will grow as we
   // discover new branch points for the cluster. The FIFO behavior is
   // accomplished by using a while-loop rather than a for-loop.  In
   // NeighborPts, the points are represented by their index in the original
   // dataset.
   // --------------------------------------------------------------------------
   var i = 0;
   while i < NeighborPts.size do {  
      // -----------------------------------------------------------------------
      // If nP was labelled NOISE during the seed search, then we know it's not
      // a branch point (it doesn't have enough neighbors), so make it a leaf
      // point of cluster C and move on.
      // -----------------------------------------------------------------------
      nP = NeighborPts[i];
      if ptlabel[nP] == -1 then {
         ptlabel[nP] = C;
      }
      // -----------------------------------------------------------------------
      // Otherwise, if nP isn't already claimed, claim it as part of C.
      // -----------------------------------------------------------------------
      else if ptlabels[nP] == 0 then {
         ptlabel[nP] = C;                         // Add nP to cluster C.
         nPNeighborPts = regionQuery(D, nP, eps); // Find all the neighbors of
                                                  // nP.
         // --------------------------------------------------------------------
         // If nP has at least MinPts neighbors, it's a branch point!  Add all
         // of its neighbors to NeighborPts by growing its domain and including
         // NpNeighborPts.
         // --------------------------------------------------------------------
         var npt = NeighborPts.size;
         var npn = NpNeighborPts.size;
         if npn  >= MinPts then {
            var nsz = npt + npn;
            dpts = {0..nsz-1};
            NeighborPts[npt..nsz-1] = NpNeightPts;
         }
      }
      i += 1;
   }
   return;
}



