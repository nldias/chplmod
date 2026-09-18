config const describe = false;
const doc = "\
================================================================================\
==> ModisTiles: simple procedures to go from the sphere to the plane            \
(and vice-versa) in a sinusoidal projection. Also identifies to which Modis tile\
a given point belongs                                                           \
================================================================================\
";
if describe then {
   writeln(doc);
   exit(0);
}
// -----------------------------------------------------------------------------
// 2026-01-23T16:14:48 a new star is born
// -----------------------------------------------------------------------------
use angles;
use ssr only amin, amax;
use Math only cos, pi;
private const RE = 6_371_007.181;  // radius of the Earth in meters
private const xmax = pi*RE;        // x maximum of plane projection
private const xmin = -xmax;        // x minimum of plane projection
private const ymax = xmax/2;       // y maximum of plane projection
private const ymin = -ymax;        // y minimum of plane projection
const DeltaMod = xmax/18;          // side of a MODIS tile
const deltapix = DeltaMod/1200;    // size of a MODIS pixel
private var
   xl,                             // x minimum of MODIS tile
   yu,                             // y maximum of MODIS tile
   xr,                             // x maximum of MODIS tile
   yl:                             // y minimum of MODIS tile
   real = nan;                     // if all are nans, ModisTileUpperLeft and/or
                                   // ModisTileLowerRight have not been called
// -----------------------------------------------------------------------------
// A relatively tight bounding box around a watershed
// -----------------------------------------------------------------------------
record RBbox {
var
   bimin,      // minimum i
   bimax,      // maximum i
   bjmin,      // minimum j
   bjmax:      // maximum j
   int;
var
   bxmin,      // minimum x of box
   bxmax,      // maximim x of box
   bymin,      // minimum y of box
   bymax:      // maximum y of box 
   real;
}
// ------------------------------------------------------------------------------
// --> FindBbox: given the boundaries, wraps a tight rectangular bounding box
// (in Sinusoidal coordinates) around  
// ------------------------------------------------------------------------------
proc FindBbox(                     // begin FindBbox                             @\label{lin:modistiles-findbbox}@ 
   const ref xb: [] real,          // The abscissas of the watershed boundary.
   const ref yb: [] real           // The ordinates of the watershed boundary.
   ): RBbox                        // A tight bounding box.
   where (xb.rank == 1 && yb.rank == 1) {                 
   const n = xb.size;
   assert ( yb.size == n ) ;
   // --------------------------------------------------------------------------
   // Maxima and minima of lons and lats.
   // --------------------------------------------------------------------------
   const maxx = amax(xb);
   const minx = amin(xb);
   const maxy = amax(yb);
   const miny = amin(yb);
   // --------------------------------------------------------------------------
   // Start looking at the corners of the region that contains the basin, with
   // some slack.
   // --------------------------------------------------------------------------
   var (i_ll,j_ll) = ijpix_xy(minx,miny);
   var (i_ul,j_ul) = ijpix_xy(minx,maxy);
   var (i_lr,j_lr) = ijpix_xy(maxx,miny);
   var (i_ur,j_ur) = ijpix_xy(maxx,maxy);
   // --------------------------------------------------------------------------
   //  Now obtain the corresponding x and y. I hope this is right!!!
   // --------------------------------------------------------------------------
   var imin = min(i_ll,i_ul,i_lr,i_ur);
   var imax = max(i_ll,i_ul,i_lr,i_ur);
   var jmin = min(j_ll,j_ul,j_lr,j_ur);
   var jmax = max(j_ll,j_ul,j_lr,j_ur);
   // --------------------------------------------------------------------------
   // The x* and y* values below are exactly on the modis pixel borders ...
   // --------------------------------------------------------------------------
   var xmin = xl + jmin*deltapix;
   var xmax = xl + (jmax+1)*deltapix;
   var ymin = yu - (imax+1)*deltapix;
   var ymax = yu - imin*deltapix;
   // ---------------------------------------------------------------------------
   // return doing nothing (for now)
   // ---------------------------------------------------------------------------
   var bbox = new RBbox(imin,imax,jmin,jmax,xmin,xmax,ymin,ymax);
   return bbox;
}  // end FindBbox                                                               @\label{lin:modistiles-endfindbbox}@
// ------------------------------------------------------------------------------
// --> ModisTileUpperLeft: define the upper left corner (in (x,y) plane
// coordinates) of the MODIS tile; store them in the private variables xl,yu
// ------------------------------------------------------------------------------
proc ModisTileUpperLeft(
   const in xin: real,
   const in yin: real
   ) {
   (xl,yu) = (xin,yin);
}
// -----------------------------------------------------------------------------
// --> ModisTileLowerRight: define the lower right corner (in (x,y) plane
// coordinates) of the MODIS tile; store them in the private variables xr,yl
// -----------------------------------------------------------------------------
proc ModisTileLowerRight(
   const in xin: real,
   const in yin: real
   ) { 
   (xr,yl) = (xin,yin);
}
// -----------------------------------------------------------------------------
// --> ModisTileCorners
// -----------------------------------------------------------------------------
proc ModisTileCorners {
   writeln("Corners of the whole tile:");
   writef("-"*20+"\n");
   writef("(%10.2dr,%10.2dr)\n",xl,yl);
   writef("(%10.2dr,%10.2dr)\n",xr,yl);
   writef("(%10.2dr,%10.2dr)\n",xl,yu);
   writef("(%10.2dr,%10.2dr)\n",xr,yu);
   writef("-"*20+"\n");
}   
// -----------------------------------------------------------------------------
// --> xySinuCoord: from (lamb,phi) to (x,y) -- sinusoidal projection to plane.
//
// Cannot say "lambda" because it is a reserved word in Chapel.
// -----------------------------------------------------------------------------
proc xySinuCoord(
   in lamb: real,             // longitude (decimals)
   in phi: real               // latitude (decimals)
   ): (real,real) {           // coordinates in the plane
   var
      x,y: real;
   phi = dec2rad(phi);        // convert phi to radians
   lamb = dec2rad(lamb);      // convert lamb to radians
   x = RE*lamb*cos(phi);
   y = RE*phi;
   return (x,y);
}

// -----------------------------------------------------------------------------
// --> geoSinuCoord: from (x,y) to (lamb,phi) (in decimals) -- inverse
// sinusoidal projection from plane to sphere (geoid)
// -----------------------------------------------------------------------------
proc geoSinuCoord(
   const in x: real,
   const in y: real
   ): (real,real) {
   var
      lamb,phi: real;
   phi = y/RE;
   lamb = x/(RE*cos(phi));
   return (rad2dec(lamb),rad2dec(phi)); // convert to degrees dec and return
}

// -----------------------------------------------------------------------------
// --> ijTile: given x and y, returns the corresponding indices (i,j) of the
// matrix of tiles.  Note that
// 
// *j is the index corresponding to x* and
// *i is the index corresponding to y* !!!
// -----------------------------------------------------------------------------
proc ijTile(
   const in x: real,     // abscissa (m) of projection onto plane
   const in y: real      // ordinate (m) of projection onto plane
   ): (int, int) {       // tile indices
   var
      x, y: real;
   assert( xmin <= x && x < xmax);
   assert( ymin < y && y <= ymax);
   var
      j,i:int;
   j = ((xmax+x)/DeltaMod):int;
   i = ((ymax-y)/DeltaMod):int;
   assert (0 <= j && j <= 35);
   assert (0 <= i && i <= 17);
   return (i,j);
}

// ------------------------------------------------------------------------------
// --> ijpix_geo: (i,j) pixel *inside a tile* from geographical (lamb,phi)
// coordinates (decimals).
// ------------------------------------------------------------------------------
proc ijpix_geo(
   const in lamb: real,
   const in phi: real
   ): (int,int) {
   var (x,y) = xySinuCoord(lamb,phi);
// -----------------------------------------------------------------------------
// Here (i,j) refers to the Modis Pixel *within* the Tile.
// -----------------------------------------------------------------------------
   var j = ((x - xl)/deltapix):int ;
   var i = ((yu - y)/deltapix):int ;
   return (i,j);
}
// ------------------------------------------------------------------------------
// --> ijpix_xy: (i,j) pixel *inside the Tile from cartesian (x,y) coordinates.
// ------------------------------------------------------------------------------
proc ijpix_xy(
   const in x: real,
   const in y: real
   ): (int,int) {
// -----------------------------------------------------------------------------
// Here (i,j) refers to the Modis Pixel *within* the Tile.
// -----------------------------------------------------------------------------
   var j = ((x - xl)/deltapix):int ;
   var i = ((yu - y)/deltapix):int ;
   return (i,j);
}
// -----------------------------------------------------------------------------
// xypix_ij: (x,y) [low,up] coordinates of (i,j) pix in the Tile.
// -----------------------------------------------------------------------------
proc xypix_ij(
   const in i: int,
   const in j: int
   ): (real,real) {
   var x = xl + j*deltapix;             // pixel's x border
   var y = yu - i*deltapix;             // pixel's y border
   return (x,y);
}
// -----------------------------------------------------------------------------
// xycenter_ij: (x,y) [center,center] coordinates of (i,j) pix in the Tile.
// -----------------------------------------------------------------------------
proc xycenter_ij(
   const in i: int,
   const in j: int
   ): (real,real) {
   var x = xl + (j+0.5)*deltapix;       // pixel's x border
   var y = yu - (i+0.5)*deltapix;       // pixel's y border
   return (x,y);
}
